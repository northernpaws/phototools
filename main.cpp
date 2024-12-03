#include <iostream>
#include <filesystem>
#include <memory>
#include <format>

#include <float.h> // FLT_MAX

#include <libraw/libraw.h>

#include <turbojpeg.h> // NOTE: is libjpeg-turbo, not TurboJPEG
#include <png.h> // libpng

// Rawspeed headers
#include <rawspeedconfig.h>
#include <io/FileReader.h>
#include <io/Buffer.h>
#include <io/FileIOException.h>
#include <parsers/RawParser.h>
#include <decoders/RawDecoder.h> // see also RafDecoder.h
#include <metadata/CameraMetadataException.h>
#include <metadata/CameraMetaData.h>

// Color management engine
#include <lcms2.h>


#include "image.h"
#include "raw_image.h"

// define this function, it is only declared in rawspeed:
#ifdef _OPENMP
// TODO: should we be using all possible cores?
extern "C" int rawspeed_get_number_of_processor_cores() { return omp_get_num_procs(); }
#else
extern "C" int __attribute__((const)) rawspeed_get_number_of_processor_cores() {
  return 1;
}
#endif

#include "src/engine/core/math/matrix.h"
#include "src/engine/core/color.h"

/* Taken from https://www.cybercom.net/~dcoffin/dcraw/dcraw.c */
/* Which itself attributes this algorithm to "Frank Markesteijn" */
const double xyz_rgb[3][3] = {			/* XYZ from RGB */
        { 0.412453, 0.357580, 0.180423 },
        { 0.212671, 0.715160, 0.072169 },
        { 0.019334, 0.119193, 0.950227 } };

const float d65_white[3] = { 0.950456, 1, 1.088754 };

#define FORC(cnt) for (c=0; c < cnt; c++)
#define FORC3 FORC(3)
#define FORC4 FORC(4)
#define FORCC FORC(colors)

#define SQR(x) ((x)*(x))
#define ABS(x) (((int)(x) ^ ((int)(x) >> 31)) - ((int)(x) >> 31))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define LIM(x,min,max) MAX(min,MIN(x,max))
#define ULIM(x,y,z) ((y) < (z) ? LIM(x,y,z) : LIM(x,z,y))
#define CLIP(x) LIM((int)(x),0,65535)
#define CLIPF(x) LIM((float)(x),0.0f,1.0f)
#define SWAP(a,b) { a=a+b; b=a-b; a=a-b; }

unsigned colors = 3; // TODO: move

// TODO: attribution to dcraw
// Convert to perceptual colorspace and differentiate in all directions.
//
// Original dcraw apparently converts to CIELab colorspace from camera, presumable from original AHD(?).
//
// Can be skipped in demosaicing to leave the image in the source colorspace.
// TODO: look at YPbPr re: darktable - https://github.com/darktable-org/darktable/blob/06f7ef340fb3010b639bbe8f85e2838a1339974e/src/iop/demosaicing/xtrans.c#L418
void cielab (const RawImage& image, ushort rgb[3], short lab[3]) {
    int c, i, j, k;

    float r, xyz[3];

    static float cbrt[0x10000], xyz_cam[3][4];

    if (!rgb) {
        for (i=0; i < 0x10000; i++) {
            r = static_cast<float>(i) / 65535.0f;
            cbrt[i] = r > 0.008856f ? pow(r,1.0f / 3.0f) : 7.787f * r + 16.0f / 116.0f;
        }

        for (i = 0; i < 3; i++)
            for (j = 0; j < colors; j++)
                for (xyz_cam[i][j] = k = 0; k < 3; k++)
                    xyz_cam[i][j] += static_cast<float>(xyz_rgb[i][k]) * image.xyz_to_cam[k][j] / d65_white[i];
        return;
    }

    xyz[0] = xyz[1] = xyz[2] = 0.5;

    FORCC {
        xyz[0] += xyz_cam[0][c] * rgb[c];
        xyz[1] += xyz_cam[1][c] * rgb[c];
        xyz[2] += xyz_cam[2][c] * rgb[c];
    }

    xyz[0] = cbrt[CLIP((int) xyz[0])];
    xyz[1] = cbrt[CLIP((int) xyz[1])];
    xyz[2] = cbrt[CLIP((int) xyz[2])];
    lab[0] = 64 * (116 * xyz[1] - 16);
    lab[1] = 64 * 500 * (xyz[0] - xyz[1]);
    lab[2] = 64 * 200 * (xyz[1] - xyz[2]);
}

#define FC(row,col) \
	(filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

/*int xtrans_fcol (const RawImage& raw_image, int row, int col) {
    return xtrans[(row+6) % 6][(col+6) % 6]
}*/

/**
 *
 * @param raw_image
 * @param row y position to get the CFA color for.
 * @param col x position to get the CFA color for.
 * @return
 */
int xtrans_fcol (const RawImage& raw_image, int row, int col) {
    // TODO: these may need to be swapped.
    uint8_t x = (col+6) % 6;
    uint8_t y = (row+6) % 6;

    const ColorFilterArray& cfa = raw_image.get_cfa();
    auto idx = y * cfa.width + x;

    return cfa.layout[idx];
}

#define TS 512		/* Tile Size */

/**
 *
 * @param raw_image The RAW image to demosaic and interpolate.
 * @param passes How many passes of the demosaic algo to run. Seems to best support 1 to 3.
 * @return The demosaiced image.
 */
// TODO: attribution to dcraw
Image xtrans_interpolate (const RawImage& raw_image, int passes) {

    unsigned short width = raw_image.get_width();
    unsigned short height = raw_image.get_height();

    if (width < TS || height < TS) {
        throw std::runtime_error("image width or height cannot be smaller then tilesize (TS) (usually 512)");
    }

    // ====================
    // Initial data conversion.
    //
    // Convert the 1-channel raw data to a 4-channel intermediate array, with the 1-channel data being distributed into
    // the first 3 channels based on its component (red, green, blue) in the color filter array for processing.

    // TODO: Even if we don't convert the algo to work on floating data, look at changing this
    //  into a more.. sane datatype, ideally with memory management (i.e. a std::vector).
    ushort (*image_data)[4];
    image_data = (ushort (*)[4]) calloc (height, width*sizeof *image_data);
    memset(image_data, 0, height * width * sizeof(*image_data));

    // TODO: also handle float
    auto raw_data = raw_image.uint16();

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // NOTE: ushort and uint16_t are the same size (2 bytes)
            auto raw_val = raw_data[y * width + x];
            auto x_fcol = xtrans_fcol(raw_image, y, x);

            image_data[y * width + x][x_fcol] = raw_val;
        }
    }

    // ===============
    // Frank Markesteijn's magic code.
    //
    // From DCRaw, modified to work with our data structures.
    // see: https://en.wikipedia.org/wiki/Dcraw

    // TODO: modify to work on float data directly to remove the intermediate conversion steps.

    int c, d, f, g, h, i, v, ng, row, col, top, left, mrow, mcol;

    int val, ndir, pass, hm[8], avg[4], color[3][8];

    static const short orth[12] = { 1,0,0,1,-1,0,0,-1,1,0,0,1 },
            patt[2][16] = { { 0,1,0,-1,2,0,-1,0,1,1,1,-1,0,0,0,0 },
                            { 0,1,0,-2,1,0,-2,0,1,1,-2,-2,1,-1,-1,1 } },
            dir[4] = { 1,TS,TS+1,TS-1 };

    short allhex[3][3][2][8], *hex;
    ushort min, max, sgrow, sgcol;
    ushort (*rgb)[TS][TS][3], (*rix)[3], (*pix)[4];
    short (*lab)    [TS][3], (*lix)[3];
    float (*drv)[TS][TS], diff[6], tr;
    char (*homo)[TS][TS], *buffer;

    fprintf (stdout, "%d-pass X-Trans interpolation...\n", passes);

//    fprintf (stdout, "X-Trans cielab...\n");
//    cielab (raw_image, 0,0);

    ndir = 4 << (passes > 1);
    buffer = (char *) malloc (TS*TS*(ndir*11+6));
    //    merror (buffer, "xtrans_interpolate()");
    rgb  = (ushort(*)[TS][TS][3]) buffer;
    lab  = (short (*)    [TS][3])(buffer + TS*TS*(ndir*6));
    drv  = (float (*)[TS][TS])   (buffer + TS*TS*(ndir*6+6));
    homo = (char  (*)[TS][TS])   (buffer + TS*TS*(ndir*10+6));

    //    fprintf (stdout, " - Map a green hexagon around each non-green pixel and vice versa\n");

    // Map a green hexagon around each non-green pixel and vice versa:
    for (row=0; row < 3; row++)
        for (col=0; col < 3; col++)
            for (ng=d=0; d < 10; d+=2) {
                g = xtrans_fcol(raw_image, row,col) == 1;
                if (xtrans_fcol(raw_image, row+orth[d],col+orth[d+2]) == 1) ng=0; else ng++;
                if (ng == 4) { sgrow = row; sgcol = col; }
                if (ng == g+1) FORC(8) {
                    v = orth[d  ]*patt[g][c*2] + orth[d+1]*patt[g][c*2+1];
                    h = orth[d+2]*patt[g][c*2] + orth[d+3]*patt[g][c*2+1];
                    allhex[row][col][0][c^(g*2 & d)] = h + v*width;
                    allhex[row][col][1][c^(g*2 & d)] = h + v*TS;
                }
            }

    //    fprintf (stdout, " - Set green1 and green3 to the minimum and maximum allowed values\n");
    // Set green1 and green3 to the minimum and maximum allowed values:
    for (row=2; row < height-2; row++)
        for (min=~(max=0), col=2; col < width-2; col++) {
            if (xtrans_fcol(raw_image, row,col) == 1 && (min=~(max=0))) continue;
            pix = image_data + row*width + col;
            hex = allhex[row % 3][col % 3][0];
            if (!max) FORC(6) {
                val = pix[hex[c]][1];
                if (min > val) min = val;
                if (max < val) max = val;
            }
            pix[0][1] = min;
            pix[0][3] = max;
            switch ((row-sgrow) % 3) {
                case 1: if (row < height-3) { row++; col--; } break;
                case 2: if ((min=~(max=0)) && (col+=2) < width-3 && row > 2) row--;
            }
        }

    for (top=3; top < height-19; top += TS-16)
        for (left=3; left < width-19; left += TS-16) {
            mrow = MIN (top+TS, height-3);
            mcol = MIN (left+TS, width-3);
            for (row=top; row < mrow; row++)
                for (col=left; col < mcol; col++)
                    memcpy (rgb[0][row-top][col-left], image_data[row*width+col], 6);
            FORC3 memcpy (rgb[c+1], rgb[0], sizeof *rgb);

            // Interpolate green horizontally, vertically, and along both diagonals:
            for (row=top; row < mrow; row++)
                for (col=left; col < mcol; col++) {
                    if ((f = xtrans_fcol(raw_image, row,col)) == 1) continue;
                    pix = image_data + row*width + col;
                    hex = allhex[row % 3][col % 3][0];
                    color[1][0] = 174 * (pix[  hex[1]][1] + pix[  hex[0]][1]) -
                                  46 * (pix[2*hex[1]][1] + pix[2*hex[0]][1]);
                    color[1][1] = 223 *  pix[  hex[3]][1] + pix[  hex[2]][1] * 33 +
                                  92 * (pix[      0 ][f] - pix[ -hex[2]][f]);
                    FORC(2) color[1][2+c] =
                                    164 * pix[hex[4+c]][1] + 92 * pix[-2*hex[4+c]][1] + 33 *
                                                                                        (2*pix[0][f] - pix[3*hex[4+c]][f] - pix[-3*hex[4+c]][f]);
                    FORC4 rgb[c^!((row-sgrow) % 3)][row-top][col-left][1] =
                            LIM(color[1][c] >> 8,pix[0][1],pix[0][3]);
                }

            for (pass=0; pass < passes; pass++) {
                if (pass == 1)
                    memcpy (rgb+=4, buffer, 4*sizeof *rgb);

                // Recalculate green from interpolated values of closer pixels:
                if (pass) {
                    for (row=top+2; row < mrow-2; row++)
                        for (col=left+2; col < mcol-2; col++) {
                            if ((f = xtrans_fcol(raw_image, row,col)) == 1) continue;
                            pix = image_data + row*width + col;
                            hex = allhex[row % 3][col % 3][1];
                            for (d=3; d < 6; d++) {
                                rix = &rgb[(d-2)^!((row-sgrow) % 3)][row-top][col-left];
                                val = rix[-2*hex[d]][1] + 2*rix[hex[d]][1]
                                      - rix[-2*hex[d]][f] - 2*rix[hex[d]][f] + 3*rix[0][f];
                                rix[0][1] = LIM(val/3,pix[0][1],pix[0][3]);
                            }
                        }
                }

                // Interpolate red and blue values for solitary green pixels:
                for (row=(top-sgrow+4)/3*3+sgrow; row < mrow-2; row+=3)
                    for (col=(left-sgcol+4)/3*3+sgcol; col < mcol-2; col+=3) {
                        rix = &rgb[0][row-top][col-left];
                        h = xtrans_fcol(raw_image, row,col+1);
                        memset (diff, 0, sizeof diff);
                        for (i=1, d=0; d < 6; d++, i^=TS^1, h^=2) {
                            for (c=0; c < 2; c++, h^=2) {
                                g = 2*rix[0][1] - rix[i<<c][1] - rix[-i<<c][1];
                                color[h][d] = g + rix[i<<c][h] + rix[-i<<c][h];
                                if (d > 1)
                                    diff[d] += SQR (rix[i<<c][1] - rix[-i<<c][1]
                                                    - rix[i<<c][h] + rix[-i<<c][h]) + SQR(g);
                            }
                            if (d > 1 && (d & 1))
                                if (diff[d-1] < diff[d])
                                    FORC(2) color[c*2][d] = color[c*2][d-1];
                            if (d < 2 || (d & 1)) {
                                FORC(2) rix[0][c*2] = CLIP(color[c*2][d]/2);
                                rix += TS*TS;
                            }
                        }
                    }

                // Interpolate red for blue pixels and vice versa:
                for (row=top+3; row < mrow-3; row++)
                    for (col=left+3; col < mcol-3; col++) {
                        if ((f = 2-xtrans_fcol(raw_image, row,col)) == 1) continue;
                        rix = &rgb[0][row-top][col-left];
                        c = (row-sgrow) % 3 ? TS:1;
                        h = 3 * (c ^ TS ^ 1);
                        for (d=0; d < 4; d++, rix += TS*TS) {
                            i = d > 1 || ((d ^ c) & 1) ||
                                ((ABS(rix[0][1]-rix[c][1])+ABS(rix[0][1]-rix[-c][1])) <
                                 2*(ABS(rix[0][1]-rix[h][1])+ABS(rix[0][1]-rix[-h][1]))) ? c:h;
                            rix[0][f] = CLIP((rix[i][f] + rix[-i][f] +
                                              2*rix[0][1] - rix[i][1] - rix[-i][1])/2);
                        }
                    }

                // Fill in red and blue for 2x2 blocks of green:
                for (row=top+2; row < mrow-2; row++) if ((row-sgrow) % 3)
                        for (col=left+2; col < mcol-2; col++) if ((col-sgcol) % 3) {
                                rix = &rgb[0][row-top][col-left];
                                hex = allhex[row % 3][col % 3][1];
                                for (d=0; d < ndir; d+=2, rix += TS*TS)
                                    if (hex[d] + hex[d+1]) {
                                        g = 3*rix[0][1] - 2*rix[hex[d]][1] - rix[hex[d+1]][1];
                                        for (c=0; c < 4; c+=2) rix[0][c] =
                                                                       CLIP((g + 2*rix[hex[d]][c] + rix[hex[d+1]][c])/3);
                                    } else {
                                        g = 2*rix[0][1] - rix[hex[d]][1] - rix[hex[d+1]][1];
                                        for (c=0; c < 4; c+=2) rix[0][c] =
                                                                       CLIP((g + rix[hex[d]][c] + rix[hex[d+1]][c])/2);
                                    }
                            }
            }
            rgb = (ushort(*)[TS][TS][3]) buffer;
            mrow -= top;
            mcol -= left;

            // Convert from camera colorspace to CIELab and differentiate in all directions:
            //
            // This converts the results to perpetual colorspace (CIELab) https://programmingdesignsystems.com/color/perceptually-uniform-color-spaces/
            //  for further processing via code.
            //
            // Perpetually uniform color spaces give a uniform appearance of brightness for the same luminance across
            // colors. This is in contrast to non-perpetual color spaces such as HSL where yellow or green tend to
            // appear brighter than other colors, even with the same lumience value.
            // see: https://programmingdesignsystems.com/color/perceptually-uniform-color-spaces/
            // TODO: this CIELab "conversion" seems to have no effect??
            //   no - putting 2 in afinity photo with a "difference" blend mode shows very, very subtle differences
            /*for (d=0; d < ndir; d++) {
                for (row=2; row < mrow-2; row++)
                    for (col=2; col < mcol-2; col++)
                        cielab (raw_image, rgb[d][row][col], lab[row][col]);
                for (f=dir[d & 3],row=3; row < mrow-3; row++)
                    for (col=3; col < mcol-3; col++) {
                        lix = &lab[row][col];
                        g = 2*lix[0][0] - lix[f][0] - lix[-f][0];
                        drv[d][row][col] = SQR(g)
                                           + SQR((2*lix[0][1] - lix[f][1] - lix[-f][1] + g*500/232))
                                           + SQR((2*lix[0][2] - lix[f][2] - lix[-f][2] - g*500/580));
                    }
            }*/

            // Build homogeneity maps from the derivatives:
            memset(homo, 0, ndir*TS*TS);
            for (row=4; row < mrow-4; row++)
                for (col=4; col < mcol-4; col++) {
                    for (tr=FLT_MAX, d=0; d < ndir; d++)
                        if (tr > drv[d][row][col])
                            tr = drv[d][row][col];
                    tr *= 8;
                    for (d=0; d < ndir; d++)
                        for (v=-1; v <= 1; v++)
                            for (h=-1; h <= 1; h++)
                                if (drv[d][row+v][col+h] <= tr)
                                    homo[d][row][col]++;
                }

            // Average the most homogenous pixels for the final result:
            if (height-top < TS+4) mrow = height-top+2;
            if (width-left < TS+4) mcol = width-left+2;
            for (row = MIN(top,8); row < mrow-8; row++)
                for (col = MIN(left,8); col < mcol-8; col++) {
                    for (d=0; d < ndir; d++)
                        for (hm[d]=0, v=-2; v <= 2; v++)
                            for (h=-2; h <= 2; h++)
                                hm[d] += homo[d][row+v][col+h];
                    for (d=0; d < ndir-4; d++)
                        if (hm[d] < hm[d+4]) hm[d  ] = 0; else
                        if (hm[d] > hm[d+4]) hm[d+4] = 0;
                    for (max=hm[0],d=1; d < ndir; d++)
                        if (max < hm[d]) max = hm[d];
                    max -= max >> 3;
                    memset (avg, 0, sizeof avg);
                    for (d=0; d < ndir; d++)
                        if (hm[d] >= max) {
                            FORC3 avg[c] += rgb[d][row][col][c];
                            avg[3]++;
                        }
                    FORC3 image_data[(row+top)*width+col+left][c] = avg[c]/avg[3];
                }
        }
    free(buffer);

    // ====================
    // Border interpolation.

    // TODO: credit dcraw
    // TODO: possibly move into it's own function if it can be re-used.
   {
        auto border = 8;
        unsigned row, col, y, x, f, c, sum[8];

        for (row=0; row < height; row++)
            for (col=0; col < width; col++) {
                if (col==border && row >= border && row < height-border)
                    col = width-border;
                memset (sum, 0, sizeof sum);
                for (y=row-1; y != row+2; y++)
                    for (x=col-1; x != col+2; x++)
                        if (y < height && x < width) {
                            f = xtrans_fcol(raw_image, y,x); // fcol
                            sum[f] += image_data[y*width+x][f];
                            sum[f+4]++;
                        }
                f = xtrans_fcol(raw_image, row,col); // fcol
                FORCC if (c != f && sum[c+4])
                        image_data[row*width+col][c] = sum[c] / sum[c+4];
            }
    }

    // ====================
    // Conversion back to float.
    //
    // Converts the 4-channel short array back to floating data after the demosaic processing.

    fprintf (stdout, "Encoding demosaiced raw to float image...\n");

    Image image = Image(width, height, 3);

    constexpr auto ushort_max = std::numeric_limits<ushort>::max();

    for (int y = 0; y < image.get_height(); y++) {
        for (int x = 0; x < image.get_width(); x++) {
            assert((y * (width*3) + (x * 3)) < image.size());
            assert(((y * (width*3) + (x * 3)) + 1) < image.size());
            assert(((y * (width*3) + (x * 3)) + 2) < image.size());

            auto pix = image_data [y * width + x];

            // Divide by max value to convert to 0.0 - 1.0 floating range.
            float r = static_cast<float>(pix[0]) / ushort_max;
            float g = static_cast<float>(pix[1]) / ushort_max;
            float b = static_cast<float>(pix[2]) / ushort_max;

            image[(y * (width*3) + (x * 3))] = r;
            image[(y * (width*3) + (x * 3)) + 1] = g;
            image[(y * (width*3) + (x * 3)) + 2] = b;
        }
    }

    // Ensure we clean up the intermediate data array.
    delete[] image_data;

    return image;
}
#undef fcol


// splits the raw data into rgb channels in the image data, doesn't look like the source image, but opening in an image
// viewer and filtering by R/G/B channels will show the (CFA) filter sensor pattern.
Image raw_colormap (const RawImage& raw_image) {
    unsigned short width = raw_image.get_width();
    unsigned short height = raw_image.get_height();

    // TODO: throw if width or height is < TS (tile size)

    ushort (*image_data)[4];
    image_data = (ushort (*)[4]) calloc (height, width*sizeof *image_data);
    memset(image_data, 0, height * width * sizeof(*image_data));

    // TODO: also handle float
    auto raw_data = raw_image.uint16();

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            assert((y * width + x) < (width*height));

            // NOTE: ushort and uint16_t are the same size (2 bytes)
            auto raw_val = raw_data[y * width + x];
            auto x_fcol = xtrans_fcol(raw_image, y, x);
            auto idx = y * width + x;

            image_data[idx][x_fcol] = raw_val;
//            image_data[y * width + x][0] = raw_val;
        }
    }

    fprintf (stdout, "Encoding raw image colormap image to float...\n");

    Image image = Image(width, height, 3);

    constexpr auto max = std::numeric_limits<ushort>::max();

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            assert((y * (width*3) + (x * 3)) < image.size());
            assert(((y * (width*3) + (x * 3)) + 1) < image.size());
            assert(((y * (width*3) + (x * 3)) + 2) < image.size());

            auto pix = image_data[y * width + x];

            // TODO: should this be l a b for CIELAB colorspace?
            float r = static_cast<float>(pix[0]) / max;
            float g = static_cast<float>(pix[1]) / max;
            float b = static_cast<float>(pix[2]) / max;

            image[(y * (width*3) + (x * 3))] = r;
            image[(y * (width*3) + (x * 3)) + 1] = g;
            image[(y * (width*3) + (x * 3)) + 2] = b;
        }
    }

    delete[] image_data;

    return image;
}

void MyLogErrorHandler(cmsContext ContextID, cmsUInt32Number ErrorCode, const char *Text) {
    printf("%s\n", Text);
}

void apply_cam_srgb_matrix(const RawImage& raw, Image& target) {
    auto embedded_matrix = raw.xyz_to_cam;

    if(target.get_channels() != 3) {
        throw std::runtime_error("can only apply matrix to 3-channel images");
    }

    core::Matrix4x3d XYZ_to_CAM;

    XYZ_to_CAM[0][0] = embedded_matrix[0][0];
    XYZ_to_CAM[0][1] = embedded_matrix[0][1];
    XYZ_to_CAM[0][2] = embedded_matrix[0][2];
    XYZ_to_CAM[1][0] = embedded_matrix[1][0];
    XYZ_to_CAM[1][1] = embedded_matrix[1][1];
    XYZ_to_CAM[1][2] = embedded_matrix[1][2];
    XYZ_to_CAM[2][0] = embedded_matrix[2][0];
    XYZ_to_CAM[2][1] = embedded_matrix[2][1];
    XYZ_to_CAM[2][2] = embedded_matrix[2][2];

    std::cout << "XYZ_to_CAM" << std::endl;
    std::cout << XYZ_to_CAM[0][0] << " " << XYZ_to_CAM[0][1] << " " << XYZ_to_CAM[0][2] << std::endl;
    std::cout << XYZ_to_CAM[1][0] << " " << XYZ_to_CAM[1][1] << " " << XYZ_to_CAM[1][2] << std::endl;
    std::cout << XYZ_to_CAM[2][0] << " " << XYZ_to_CAM[2][1] << " " << XYZ_to_CAM[2][2] << std::endl;

    // sRGB D65?
    static const core::Matrix3x3d RGB_to_XYZ {
            { 0.412453, 0.357580, 0.180423 },
            { 0.212671, 0.715160, 0.072169 },
            { 0.019334, 0.119193, 0.950227 },
    };

    std::cout << "RGB_to_XYZ" << std::endl;
    std::cout << RGB_to_XYZ[0][0] << " " << RGB_to_XYZ[0][1] << " " << RGB_to_XYZ[0][2] << std::endl;
    std::cout << RGB_to_XYZ[1][0] << " " << RGB_to_XYZ[1][1] << " " << RGB_to_XYZ[1][2] << std::endl;
    std::cout << RGB_to_XYZ[2][0] << " " << RGB_to_XYZ[2][1] << " " << RGB_to_XYZ[2][2] << std::endl;

    // adapted from https://github.com/darktable-org/darktable/blob/c6e02200a15bee75557f4f98a8cbbf27c8b2157f/src/common/colorspaces.c#L2317
    core::Matrix4x3d RGB_to_CAM;

    // Multiply RGB matrix
    for(int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            RGB_to_CAM[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
                RGB_to_CAM[i][j] += XYZ_to_CAM[i][k] * RGB_to_XYZ[k][j];
        }
    }

    // Normalize cam_rgb so that cam_rgb * (1,1,1) is (1,1,1,1)
    for(int i = 0; i < 4; i++) {
        double num = 0.0;
        for(int j = 0; j < 3; j++)
            num += RGB_to_CAM[i][j];
        for(int j = 0; j < 3; j++)
            RGB_to_CAM[i][j] /= num;
        // TODO: important? may be able to provide camera whitebalance multipliers here to do matrix application and whitebalance correction in one step
//        if(mul) mul[i] = 1.0 / num;
    }

    std::cout << "RGB_to_CAM" << std::endl;
    std::cout << RGB_to_CAM[0][0] << " " << RGB_to_CAM[0][1] << " " << RGB_to_CAM[0][2] << std::endl;
    std::cout << RGB_to_CAM[1][0] << " " << RGB_to_CAM[1][1] << " " << RGB_to_CAM[1][2] << std::endl;
    std::cout << RGB_to_CAM[2][0] << " " << RGB_to_CAM[2][1] << " " << RGB_to_CAM[2][2] << std::endl;

    core::Matrix3x3d out_CAM_to_RGB;

    // TODO: add method for extracting 3x3 from 4x3
    out_CAM_to_RGB = core::inverse<double>({
         {RGB_to_CAM[0][0], RGB_to_CAM[0][1], RGB_to_CAM[0][2]},
         {RGB_to_CAM[1][0], RGB_to_CAM[1][1], RGB_to_CAM[1][2]},
         {RGB_to_CAM[2][0], RGB_to_CAM[2][1], RGB_to_CAM[2][2]}
    });

    std::cout << "out_CAM_to_RGB" << std::endl;
    std::cout << out_CAM_to_RGB[0][0] << " " << out_CAM_to_RGB[0][1] << " " << out_CAM_to_RGB[0][2] << std::endl;
    std::cout << out_CAM_to_RGB[1][0] << " " << out_CAM_to_RGB[1][1] << " " << out_CAM_to_RGB[1][2] << std::endl;
    std::cout << out_CAM_to_RGB[2][0] << " " << out_CAM_to_RGB[2][1] << " " << out_CAM_to_RGB[2][2] << std::endl;

    target.color_transform(out_CAM_to_RGB);
}

Image to_srgb(Image& source) {
    cmsSetLogErrorHandler(MyLogErrorHandler);

    cmsHPROFILE hInProfile, hOutProfile; cmsHTRANSFORM hTransform;

    // Load the applicable colorspaces
//    hInProfile = cmsOpenProfileFromFile("Generic Lab Profile.icc", "r"); // Generic Lab Profile //Lab-D50_2deg.icc
    hInProfile = cmsCreateLab4Profile(cmsD50_xyY());

//    hInProfile = cmsCreateXYZProfile();

    if (hInProfile == nullptr) {
        throw std::runtime_error("Failed to load CIELab color profile!");
    }

//    hOutProfile = cmsOpenProfileFromFile("sRGB Profile.icc", "r");
    hOutProfile = cmsCreate_sRGBProfile();
    if (hOutProfile == nullptr) {
        throw std::runtime_error("Failed to load sRB color profile!");
    }
//cmsSigLabData
    hTransform = cmsCreateTransform(hInProfile, TYPE_Lab_FLT, // TYPE_RGB_FLT,
                                    hOutProfile, TYPE_RGB_FLT,
                                    INTENT_PERCEPTUAL,
                                    0); // flag cmsFLAGS_BLACKPOINTCOMPENSATION?
    if (hTransform == nullptr) {
        throw std::runtime_error("Failed to create color profile transform!");
    }

    cmsCloseProfile(hInProfile);
    cmsCloseProfile(hOutProfile);

    // Create a new image with the output colorspace.
    Image dst = Image(source.get_width(), source.get_height(), source.get_channels());

    float min = 0;
    float max = 0;

    for (int y = 0; y < source.get_height(); y++) {
        for (int x = 0; x < source.get_width(); x++) {
            auto val = source[y * source.get_width() + x];

            if (val < min) {
                min = val;
            }

            if (val > max) {
                max = val;
            }
        }
    }

    std::cout << "input min: " << min << std::endl;
    std::cout << "input max: " << max << std::endl;

    // Transform the source image to the dst image.
    //
    // NOTE: Size is width*height in pixels, excluding channels. Channels
    //  are compensated for based on the cms input and output format.
    cmsDoTransform(hTransform, source.get_data(), dst.get_data(), source.get_width() * source.get_height());

    cmsDeleteTransform(hTransform);

    min = 0;
    max = 0;

    for (int y = 0; y < dst.get_height(); y++) {
        for (int x = 0; x < dst.get_width(); x++) {
            auto val = dst[y * dst.get_width() + x];

            if (val < min) {
                min = val;
            }

            if (val > max) {
                max = val;
            }
        }
    }

    std::cout << "output min: " << min << std::endl;
    std::cout << "output max: " << max << std::endl;

    // TODO: do we need to correct for values <0 >1?

    return dst;
}

// gamma = 2.4 for srgb according to rawtherapee
// slope = 12.92  for srgb according to rawtherapee
//
//
// ref: toe curve https://www.dpreview.com/forums/thread/4466218
// TODO: is it called slope or toe?
void apply_gamma_curve(Image& image, float gamma, float toe) {
    // TODO: Implementations online seem to like to use a LUT for speed with the pow calculation, but that's difficult
    //  if we apply to curve when we're using floating-point data. Maybe we should apply a base gamma curve before
    //  converting the RAW data from uint16 to float?
    // Built a LUT (lookup table) for converting
//    unsigned char lut[256];
//    for (int i = 0; i < 256; i++) {
//        lut[i] = pow((float)(i / 255.0), gamma) * 255.0f);
//    }
    auto data = image.get_data();
    for (int i = 0; i < image.get_height(); i++) {
        for (int j = 0; j < image.get_width(); j++) {
            float *pix = &(data[image.pixel_index(j, i)]);
            for (int c = 0; c < image.get_channels(); c++) {
                // For gamma with slope
                // ref: https://stackoverflow.com/a/13558570
                // see: https://stackoverflow.com/a/17343790 ( PHP implementation of the WCAG 2.0 SC 1.4.3 relative luminance and contrast ratio formulae)
                if(pix[c] <= 0.0031308f) { // TODO: document what specification  these thresholds are from
                    pix[c] *= toe;
                } else {
                    // TODO: multiple sources say the following should be correct, but it doesn't look right
                    pix[c] = 1.055f * pow(pix[c], 1.0f / gamma) - 0.055f;
                }

                // clamp the low-end.
                // TODO: see if this is actually needed - all the examples using opencv used cv::saturate_cast so maybe?
                if (pix[c] < 0) {
                    pix[c] = 0;
                }

                // For basic gamma
                /*pix[c] = pow(pix[c], 1/gamma);

                // clamp the low-end.
                // TODO: see if this is actually needed - all the examples using opencv used cv::saturate_cast so maybe?
                if (pix[c] < 0) {
                    pix[c] = 0;
                }*/
            }
        }
    }
}

void apply_whitebalance(Image& image, std::array<float, 4> wb) {
    auto channels = image.get_channels();
    auto stride = image.get_width() * channels;

    for (int y = 0; y < image.get_height(); ++y) {
        for (int x = 0; x < image.get_width(); ++x) {
            auto pix = (y * stride + (x * channels));
            assert(pix < (image.get_width() * image.get_height() * image.get_channels()));

            // TODO: also apply the 4th value if needed
            for (int c = 0; c < 3; ++c) {
                auto idx = pix + c;

                assert(idx < (image.get_width() * image.get_height() * image.get_channels()));

                image[idx] = image[idx] * wb[c];
            }
        }
    }
}

float decimal_value(int numerator, int denominator) {
    float decimalValue = 0;

    if (denominator > 0) {
        decimalValue = static_cast<float>(numerator) / static_cast<float>(denominator);
    }

    return decimalValue;
}

RawImage decode_image_rawspeed(const char* filename) {
    std::unique_ptr<rawspeed::CameraMetaData> cameras = nullptr;

    // #ifdef HAVE_PUGIXML
    // Initialize the camera metadata bank.
    if (cameras == nullptr) {
        std::cout << "loading camera metadata" << std::endl;
        try {
            cameras = std::make_unique<rawspeed::CameraMetaData>("cameras.xml");
        } catch (rawspeed::CameraMetadataException &e) {
            // Reading metadata failed. e.what() will contain error message.
            std::cout << "ERROR failed to load camera metadata: " << e.what() << std::endl;
        }
    }
//#endif

    rawspeed::FileReader reader(filename);

    std::unique_ptr<std::vector<
            uint8_t, rawspeed::DefaultInitAllocatorAdaptor<
                    uint8_t, rawspeed::AlignedAllocator<uint8_t, 16>>>>
            storage;
    rawspeed::Buffer buffer;

    try {
        std::tie(storage, buffer) = reader.readFile();
    } catch (rawspeed::FileIOException &e) {
        // Handle errors
        std::cout << "raw file io exception " << e.what() << std::endl;
        return RawImage{};
    }

    // Parse the file data buffer.
    rawspeed::RawParser parser(buffer);

    // Create a decoder derived from the parser.
    std::unique_ptr<rawspeed::RawDecoder> decoder = parser.getDecoder(cameras.get());
    if (decoder == nullptr) {
        std::cout << "couldn't get decoder " << std::endl;
        return RawImage{};
    }

    decoder->failOnUnknown = true; // refuses to decode unknown camera, even if the data can still be decoded
    decoder->applyCrop = true; // crop the raw to the active image area (as opposed to the entire sensor area)

    try {
        // throws RawDecoderException on issue
        decoder->checkSupport(cameras.get()); // check if the decoder supports the camera in the metadata
    } catch(rawspeed::RawDecoderException &e) {
        // Handle decoder error
        std::cout << "raw decoder exception " << e.what() << std::endl;
        return RawImage{};
    }

    // NOTE: You will most likely find that a relatively long time is spent actually reading the file. The
    //   biggest trick to speeding up raw reading is to have some sort of prefetching going on while the
    //   file is being decoded. This is the main reason why RawSpeed decodes from memory, and doesn’t use
    //   direct file reads while decoding.
    //
    // The simplest solution is to start a thread that simply reads the file, and rely on the system cache
    //   to cache the data. This is fairly simple and works in 99% of all cases. So if you are doing batch
    //   processing simply start a process reading the next file, when the current image starts decoding.
    //   This will ensure that your file is read linearly, which gives the highest possible throughput.

    // TODO: manual bad-pixel interpolation: https://github.com/darktable-org/rawspeed/tree/develop/src/librawspeed#bad-pixel-elimination

    // Decode the image
    try {
        decoder->decodeRaw();
    } catch(rawspeed::RawDecoderException &e) {
        // Handle decoder error
        std::cout << "raw decoder exception " << e.what() << std::endl;
    }

    // Retrieves metadata from the cameras.xml file to apply adjustments such as sensor
    // cropping, black-point adjustment, and default CFA pattern to the decoded raw.
    decoder->decodeMetaData(cameras.get());

    // NOTE: following functions operate on the image data after it's
    // been cropped to the active image area from the sensor area.

    // untouched raw data, cropped to the active image area from metadata
    // throws RawDecoderException on error
    rawspeed::RawImage raw = decoder->mRaw;
    // TODO: do we need to check if raw is undefined (see mRaw comment)?

    // Decode for decoding errors.
    for (auto const& error: raw->getErrors()) {
        std::cout << "RAWSPEED ERROR: " << error << std::endl;
    }

    // decoder->fujiRotate

//    if(!raw->blackAreas.empty()) {
//        raw->calculateBlackAreas();
//    }

    // TODO: remove and do scaling separately
    // apply black/white scaling to the image, data is normalized into the 0->65535 range for 16-bit images
    raw->scaleBlackWhite();

    // get raw image data

    // 1 on CFA images, and 3, found in some DNG images.
    // cannot assume that an images is CFA just because it is 1 cpp - greyscale dng images from things like scanners can be saved like that
    uint32_t components_per_pixel = raw->getCpp();
    std::cout << "components_per_pixel " << components_per_pixel << "(cpp)" << std::endl;

    if (components_per_pixel == 1) {
        std::cout << "image is probably a CFA or grayscale scan" << std::endl;
    }

    // TYPE_USHORT16 (most common) - unsigned 16 bit data
    // TYPE_FLOAT32 (found in some DNGs)
    rawspeed::RawImageType type = raw->getDataType();
    if (type == rawspeed::RawImageType::UINT16) {
        std::cout << "type: TYPE_USHORT16" << std::endl;
    } else if (type == rawspeed::RawImageType::F32) {
        std::cout << "type: TYPE_FLOAT32" << std::endl;
    }

    // indicates whether the image has all components per pixel, or if it was taken with a colorfilter array
    // usually corresponds to the number of components per pixel (1 on CFA, 3 on non-CFA).
    bool is_cfa = raw->isCFA;


    // image data
    // ref: https://github.com/darktable-org/rawspeed/issues/491
    // size of array should be width * height * bytes_per_pixel
    auto data2d = raw->getU16DataAsCroppedArray2DRef();
//    uint16_t* data = &(data2d(0, 0)); //  pointer to the Raw data at pixel x, y, coordinate after crop, Don’t use this function on every pixel, but instead increment the pointer yourself.

    // size should be
    int width = raw->dim.x;
    int height = raw->dim.y;
    int pitch_in_bytes = raw->pitch; // number of bytes between lines, since this is usually NOT width * components_per_pixel * bytes_per_component. calculate a pointer at line y, use &data[y * raw->pitch].

    uint32_t bytes_per_pixel = raw->getBpp();

    std::cout << "image width " << width << std::endl;
    std::cout << "image height " << height << std::endl;

    std::cout << "image pitch " << pitch_in_bytes << "(bytes between lines)" << std::endl;
    std::cout << "image bytes per pixel(component?) " << bytes_per_pixel <<  std::endl;

    const rawspeed::Camera* cam = cameras->getCamera(raw->metadata.make, raw->metadata.model, raw->metadata.mode);
    if (cam->cropAvailable) {
        std::cout << "vendor camera crop available" << std::endl;
        std::cout << "crop pos: " << cam->cropPos.x << "x" << cam->cropPos.y << " size: " << cam->cropSize.x << "x" << cam->cropSize.y << std::endl;
    } else {
        std::cout << "vendor camera crop not available" << std::endl;
    }

    // ====================
    // Decode the Color Filter Array (CFA) information.

    ColorFilterArray raw_image_cfa;
    if (is_cfa) {
        rawspeed::ColorFilterArray cfa = raw->cfa;

        std::cout << "raw CFA:" << std::endl;

        std::cout << cfa.asString() << std::endl;

        std::cout << "correcting CFA for crop " << raw->getCropOffset().x << "x" << raw->getCropOffset().y << std::endl;

        // NOTE X-Trans CFA is not corrected, but is corrected for other sensors.
        //  see: https://github.com/darktable-org/rawspeed/issues/498#issuecomment-1624448605
        if (cfa.getSize().x == 6 && cfa.getSize().y == 6) {
            /*
             * GREEN,GREEN,RED,GREEN,GREEN,BLUE
             * GREEN,GREEN,BLUE,GREEN,GREEN,RED
             * BLUE,RED,GREEN,RED,BLUE,GREEN
             * GREEN,GREEN,BLUE,GREEN,GREEN,RED
             * GREEN,GREEN,RED,GREEN,GREEN,BLUE
             * RED,BLUE,GREEN,BLUE,RED,GREEN
             */

            // TODO: need to figure out if this is a determinable value or not for xtrans
            // TODO: is the need to apply an offset because this image is rotated? does the cfa need to be rotated???

            // was originally 1x3 or 4x0, but I think the matrix was rotated oddly then??
            // now appears to be 3x4
            cfa.shiftRight(3); // 4
            cfa.shiftDown(4);//0
        }

        std::cout << "corrected CFA:" << std::endl;

        std::cout << cfa.asString() << std::endl;

        auto size = cfa.getSize();

        int cfa_width = size.x;
        int cfa_height = size.y;

        auto crop_x = raw->getCropOffset().x;
        auto crop_y = raw->getCropOffset().y;

        raw_image_cfa.width = cfa_width;
        raw_image_cfa.height = cfa_height;

        raw_image_cfa.layout = std::vector<unsigned char>(cfa_width * cfa_height);

        // Compatability conversion to tiff CFAPattern format for the demosaic algo.
        // https://www.awaresystems.be/imaging/tiff/tifftags/privateifd/exif/cfapattern.html
        for (int y = 0; y < cfa_height; y++) {
            for (int x = 0; x < cfa_width; x++) {
                assert((y * cfa_width + x) < cfa_width * cfa_height);

                rawspeed::CFAColor c = cfa.getColorAt(x,y);

                // NOTE: rawspeed::CFAColor matches the TIFF CFAPattern specification.
//                ret_image.xtrans[x][y] = static_cast<char>(c);

                // TODO: convert to concrete enum type
                raw_image_cfa.layout[y * cfa_width + x] = static_cast<unsigned char>(c);
            }
        }
    }

    // load data and cast to float
    // TODO: opemmp parallelization

    // ===============
    // Copy the RAW image data from the Rawspeed handle to our own buffer
    // to preserve ownership over the pointer lifecycle. Otherwise, Rawspeed
    // frees the data automatically once we're finished with the decoder.

    // TODO: use vector instead
    // = new uint16_t[width*height]
    // TODO: should this be just width*height? memory issues?
    std::unique_ptr<uint16_t[]> raw_data = std::make_unique<uint16_t[]>(width*height);

    // TODO: look at using memcpy to do the copy in one go.
    // NOTE: not sure memcopy is reliable as the source data for the cropped array is the uncropped data?
    //auto src_ptr = &(raw->getU16DataAsCroppedArray2DRef()(0, 0));
    //memcpy(ret_image.raw_data, src_ptr, components_per_pixel*width*height*sizeof(uint16_t));

    for(int y = 0; y < height; y++) {
        auto ptr = data2d[y].begin();

        for (int x = 0; x < width; x++) {
            assert((y * width + x) < (width * height));

            raw_data.get()[y * width + x] = *(ptr + x);
        }
    }

    RawImage raw_image = RawImage(UINT16, width, height, raw_image_cfa, std::move(raw_data));

    // ====================
    // Decode the camera color matrix(es)

    std::cout << "colorMatrix.capacity: " << raw->metadata.colorMatrix.capacity() << std::endl;

    // TODO: load matrix into new RawImage class
    // row-major color matrix from RAW metadata to convert XYZ values to camera native color space under D65.
    auto mat = raw->metadata.colorMatrix;

    std::cout << "rational xyz_to_cam matrix size " << mat.size() << std::endl;

    // TODO: these are different from what a Adobe converted DNG outputs?
    std::cout << "rational xyz_to_cam matrix:" << std::endl;
    std::cout << std::format("{}/{} {}/{} {}/{}", mat[0].num, mat[0].den, mat[1].num, mat[1].den, mat[2].num, mat[2].den) << std::endl;
    std::cout << std::format("{}/{} {}/{} {}/{}", mat[3].num, mat[3].den, mat[4].num, mat[4].den, mat[5].num, mat[5].den) << std::endl;
    std::cout << std::format("{}/{} {}/{} {}/{}", mat[6].num, mat[6].den, mat[7].num, mat[7].den, mat[8].num, mat[8].den) << std::endl;

    int i = 0;
    for (int y = 0; y < 3; y++) {
        for (int x = 0; x < 3; x++) {
            // TODO: are x and y swapped? - yes, they were...
            raw_image.xyz_to_cam[y][x] = decimal_value(mat[i].num, mat[i].den);

            i++;
        }
    }

    std::cout << "decimal xyz_to_cam matrix:" << std::endl;
    std::cout << raw_image.xyz_to_cam[0][0] << " " << raw_image.xyz_to_cam[0][1] << " " << raw_image.xyz_to_cam[0][2] << std::endl;
    std::cout << raw_image.xyz_to_cam[1][0] << " " << raw_image.xyz_to_cam[1][1] << " " << raw_image.xyz_to_cam[1][2] << std::endl;
    std::cout << raw_image.xyz_to_cam[2][0] << " " << raw_image.xyz_to_cam[2][1] << " " << raw_image.xyz_to_cam[2][2] << std::endl;

    // monochrome matrix for debugging: https://github.com/Beep6581/RawTherapee/blob/c72e67ae24b00cb669165a965cbd6633ff44e4de/rtengine/rawimage.cc#L358
//    ret_image.rgb_cam[0][0] = 1; ret_image.rgb_cam[1][0] = 0; ret_image.rgb_cam[2][0] = 0;
//    ret_image.rgb_cam[0][1] = 0; ret_image.rgb_cam[1][1] = 1; ret_image.rgb_cam[2][1] = 0;
//    ret_image.rgb_cam[0][2] = 0; ret_image.rgb_cam[1][2] = 0; ret_image.rgb_cam[2][2] = 1;


    std::cout << "as shot whitebalance coefficients: " << raw->metadata.wbCoeffs[0] << " " << raw->metadata.wbCoeffs[1] << " " << raw->metadata.wbCoeffs[2]<< " " << raw->metadata.wbCoeffs[3] << std::endl;

    auto wb_r = raw->metadata.wbCoeffs[0];
    auto wb_g = raw->metadata.wbCoeffs[1];
    auto wb_b = raw->metadata.wbCoeffs[2];

    // Divide the R and B values by G to get the whitebalance multipliers.
    // NOTE: white is g the anchor/reference?
    // TODO: Darktable calculates a temp and tint from these as well
    raw_image.whitebalance[0] = wb_r / wb_g;
    raw_image.whitebalance[1] = wb_g / wb_g;
    raw_image.whitebalance[2] = wb_b / wb_g;


    raw_image.white_point = raw->whitePoint.value_or(65535); // use uint16 max if no whitepoint was defined
    raw_image.black_level = raw->blackLevel;

    // TODO: calculate seperate black levels: https://github.com/darktable-org/darktable/blob/c6e02200a15bee75557f4f98a8cbbf27c8b2157f/src/imageio/imageio_rawspeed.cc#L211

    std::cout << "converted whitebalance coefficient mltipliers: " << raw_image.whitebalance[0] << " " << raw_image.whitebalance[1] << " " << raw_image.whitebalance[2]<< " " << raw_image.whitebalance[3] << std::endl;

    std::cout << "decoded raw image!" << std::endl;

    return raw_image;
}


void write_png_raw(const RawImage& raw_image, const std::string& filename) {
    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) abort();

    // TODO: populate warning func pointers
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) abort();

    png_infop info = png_create_info_struct(png);
    if (!info) abort();

    // injects standard error handling if no error handlers where provided to the create_write_struct call
    if (setjmp(png_jmpbuf(png))) abort();

    // use internal streaming write methods.
    png_init_io(png, fp);

    std::cout << "writing 1 channel (grayscale) png" << std::endl;
// Output is 8bit depth, grayscale format.
    png_set_IHDR(
            png,
            info,
            raw_image.get_width(), raw_image.get_height(),
            8, // bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on the color_type selected
            PNG_COLOR_TYPE_GRAY,
            PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT
    );

    png_write_info(png, info);

    png_uint_32 height, width;

    height = raw_image.get_height();
    width = raw_image.get_width();

    int rowbytes = png_get_rowbytes(png, info);

    std::cout << "width: " << width << " height: " << height << " row bytes: " << rowbytes << std::endl;

    std::vector<png_byte> data(rowbytes * height);
    std::vector<png_bytep> row_pointers(height);

    // TODO: support float raws.
    uint16_t* raw_data = raw_image.uint16();

    for (size_t i = 0; i < height; ++i) {
        row_pointers[i] = (unsigned char *) &data[i * rowbytes];
    }
    std::cout << "writing 1 channel (grayscale) png data..." << std::endl;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            auto ptr = data.data() + (y * width + x);

            // NOTE: we assume single-channel for RAW image
            // TODO: we can change the PNG depth to 16 bit and avoid this nasty converting
            unsigned char val = static_cast<unsigned char>((static_cast<float>(raw_data[y * width + x]) / std::numeric_limits<uint16_t>::max()) * std::numeric_limits<unsigned char>::max());  // NOTE: 8 bit (255)

            *ptr = val;
        }
    }

    // write the image in one go

    // To remove the alpha channel for PNG_COLOR_TYPE_RGB format,
    // Use png_set_filler().
    //png_set_filler(png, 0, PNG_FILLER_AFTER);

    png_write_image(png, row_pointers.data());
    png_write_end(png, NULL);

    // close the file and cleanup
    fclose(fp);


    png_destroy_write_struct(&png, &info);
}

void write_png_image(const Image& image, std::string filename) {
    FILE *fp = fopen(filename.c_str(), "wb");
    if (!fp) abort();

    // TODO: populate warning func pointers
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) abort();

    png_infop info = png_create_info_struct(png);
    if (!info) abort();

    // injects standard error handling if no error handlers where provided to the create_write_struct call
    if (setjmp(png_jmpbuf(png))) abort();

    // use internal streaming write methods.
    png_init_io(png, fp);

//    std::cout << "using linear gamma" << std::endl;
    // Set to linear gamma instead of sRGB default (TODO: make option later depending on image colorspace)
//    png_set_gamma(png, PNG_GAMMA_LINEAR, PNG_GAMMA_LINEAR);

    if (image.get_channels() == 1) {
        std::cout << "writing 1 channel (grayscale) png" << std::endl;
        // Output is 8bit depth, grayscale format.
        png_set_IHDR(
                png,
                info,
                image.get_width(), image.get_height(),
                8, // bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on the color_type selected
                PNG_COLOR_TYPE_GRAY,
                PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT
        );
    } else if (image.get_channels() == 3) {
        std::cout << "writing 3 channel (RGB) png" << std::endl;
        // Output is 8bit depth, RGB format.
        png_set_IHDR(
                png,
                info,
                image.get_width(), image.get_height(),
                8, // bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on the color_type selected
                PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT
        );
    } else {
        throw std::runtime_error("unsupported channels for png export");
    }


    png_write_info(png, info);

    // TODO: embed icc profile
//    png_set_iCCP(png, info, )

    png_uint_32 height, width;
    height = image.get_height();
    width = image.get_width();

    int rowbytes = png_get_rowbytes(png, info);

    std::cout << "width: " << width << " height: " << height << " row bytes: " << rowbytes << std::endl;

    std::vector<png_byte> data(rowbytes * height);
    std::vector<png_bytep> row_pointers(height);
    for (size_t i = 0; i < height; ++i) {
        row_pointers[i] = (unsigned char *) &data[i * rowbytes];
    }

    if (image.get_channels() == 1) {
        std::cout << "writing 1 channel (grayscale) png data..." << std::endl;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                auto ptr = data.data() + (y * width + x);

                // NOTE: we assume single channel for RAW image
                auto val = static_cast<unsigned char>(image[y * width + x] * 255.0f); // NOTE: 8 bit (255)

                *ptr = val;
            }
        }
    } else if (image.get_channels() == 3) {
        std::cout << "writing 3 channel (RGB) png data..." << std::endl;
//        auto data = image.get_data();
        for(int y = 0; y < height; y++) {
            png_bytep row = row_pointers[y];
            for(int x = 0; x < width; x++) {
                png_bytep px = &(row[x * 3]); // 3 = 3 channels (RGB)
                const float* dat_px = &(image[image.pixel_index(x, y)]);

                auto r = static_cast<unsigned char>(dat_px[0] * 255.0f);  // NOTE: 8 bit (255);
                auto g = static_cast<unsigned char>(dat_px[1] * 255.0f);  // NOTE: 8 bit (255);
                auto b = static_cast<unsigned char>(dat_px[2] * 255.0f);  // NOTE: 8 bit (255);

                px[0] = r;
                px[1] = g;
                px[2] = b;
            }
        }
    } else {
        throw std::runtime_error("unsupported channels for png export");
    }


    // write the image in one go

    // To remove the alpha channel for PNG_COLOR_TYPE_RGB format,
    // Use png_set_filler().
    //png_set_filler(png, 0, PNG_FILLER_AFTER);

    png_write_image(png, row_pointers.data());
    png_write_end(png, NULL);

    // close the file and cleanup
    fclose(fp);

    png_destroy_write_struct(&png, &info);
}

int main() {
    std::cout << "Current path is " << std::filesystem::current_path() << '\n';

    // we use float32 processing in the image pipeline.
    static_assert(sizeof(float) * CHAR_BIT == 32, "require 32 bits floats");

    auto filename = "DSCF1374.RAF"; // "DSCF1785.RAF"; // "DSCF1374.RAF";

    try {
        std::cout << "decoding raw image...\n";
        auto raw_image = decode_image_rawspeed(filename);

        std::cout << "checking image validity...\n";
        if (!raw_image.is_valid()) {
            std::cout << "invalid raw image from decoder" << std::endl;
            return -1;
        }

        // Write out test image
        /*std::cout << "writing test image...\n";
        write_png_raw(raw_image, "raw_data.png");*/

        // Test image showing pixel color decodings
        /*std::cout << "creating raw colormap test image..." << std::endl;
        auto color_image = raw_colormap(raw_image);
        write_png_image(color_image, "raw_colormap.png");*/

        // Preform demosaicing on the raw image data.
        std::cout << "demosaicing..." << std::endl;
        auto image = xtrans_interpolate(raw_image, 3);

        // Write a test image showing the demosaiced raw image (with no furthur gamma or colorspace correction)
        //std::cout << "writing demosaiced image..." << std::endl;
        //write_png_image(image, "demosaiced_data.png");

        std::cout << "applying whitebalance correction to demosaiced image..." << std::endl;
        apply_whitebalance(image, raw_image.whitebalance);
        //write_png_image(image, "demosaiced_data_wb_corrected.png");

        std::cout << "applying camera srgb correction matrix..." << std::endl;
        apply_cam_srgb_matrix(raw_image, image);
        //write_png_image(image, "demosaiced_data_wb_corrected_srgb.png");

        std::cout << "applying srgb gamma curve matrix..." << std::endl;
        apply_gamma_curve(image, 2.4f, 12.92f);
        write_png_image(image, "demosaiced_data_wb_corrected_srgb_gamma.png");
    } catch (const std::exception& ex) {
        std::cout << "decoding exception " << ex.what() << std::endl;
    }

    return 0;
}
