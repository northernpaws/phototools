
## Exif tag data:
```
exiv2 -pt DSCF1374.RAF
```
Or:
```
exiftool DSCF1374.RAF
```

# DCRAW

## Extract RAW data from RAW image:
```
dcraw -4 -D -T <filename>
```

## Convert with camera "as shot" whitebalance:
```
dcraw -o 1 -w DSCF1374.RAF
```

# Colorspaces

## Extract with libraw:

```
$ brew install libraw
$ raw-identify -v DSCF1374.RAF 
...
Makernotes WB data:               coeffs                  EVs
  As shot                   346 302 912 0    0.20  0.00  1.59  0.00
  Tungsten (Incandescent)    349  302  825  302    0.21  0.00  1.45  0.00
  Fine Weather               524  302  538  302    0.80  0.00  0.83  0.00
  Shade                      572  302  460  302    0.92  0.00  0.61  0.00
  Daylight Fluorescent       685  302  466  302    1.18  0.00  0.63  0.00
  Day White Fluorescent      582  302  580  302    0.95  0.00  0.94  0.00
  Cool White Fluorescent     555  302  739  302    0.88  0.00  1.29  0.00
  Illuminant A               349  302  825  302    0.21  0.00  1.45  0.00
  D65                        598  302  471  302    0.99  0.00  0.64  0.00
  Camera Auto                346  302  912  302    0.20  0.00  1.59  0.00

Camera2RGB matrix (mode: 1):
1.2582  0.0530  -0.3112
-0.1169 1.5897  -0.4728
-0.0984 -0.5063 1.6046

XYZ->CamRGB matrix:
1.3426  -0.6334 -0.1177
-0.4244 1.2136  0.2371
0.0580  0.1303  0.5980

```

## Calculating camera colorspace matrixes:

https://github.com/wjakob/hdrmerge/blob/1dbb0d2b942343567cda721863cd86683635fd80/main.cpp#L269C12-L288C26
```
<< "*******************************************************************************" << endl
 << "Warning: no sensor2xyz matrix was specified -- this is necessary to get proper" << endl
 << "sRGB / XYZ output. To acquire this matrix, convert any one of your RAW images" << endl
 << "into a DNG file using Adobe's DNG converter on Windows / Mac (or on Linux," << endl
 << "using the 'wine' emulator). The run" << endl
 << endl
 << "  $ exiv2 -pt the_image.dng 2> /dev/null | grep ColorMatrix2" << endl
 << "  Exif.Image.ColorMatrix2 SRational 9  <sequence of ratios>" << endl
 << endl
 << "The sequence of a rational numbers is a matrix in row-major order. Compute its" << endl
 << "inverse using a tool like MATLAB or Octave and add a matching entry to the" << endl
 << "file hdrmerge.cfg (creating it if necessary), like so:" << endl
 << endl
 << "# Sensor to XYZ color space transform (Canon EOS 50D)" << endl
 << "sensor2xyz=1.933062 -0.1347 0.217175 0.880916 0.725958 -0.213945 0.089893 " << endl
 << "-0.363462 1.579612" << endl
 << endl
 << "-> Providing output in the native sensor color space, as no matrix was given." << endl
 << "*******************************************************************************" << endl
 << endl;
```
Example on test file:
```
# first use Adobe Digital Negative Converted to convert RAF file
$ cd dng_test/
$ exiv2 -pt DSCF1374.DNG 2> /dev/null | grep ColorMatrix2 
Exif.Image.ColorMatrix2                      SRational   9  13717/10000 -6490/10000 -1154/10000 -4348/10000 12266/10000 2335/10000 -690/10000 1286/10000 6134/10000
```
Gives us the rational ratios matrix for XYZ to camera sensor colorspace of:
```
 13717/10000 -6490/10000  -1154/10000
-4348/10000   12266/10000  2335/10000
-690/10000    1286/10000   6134/10000
```
In decimals:
```
1.3717   -0.649  -0.1154
-0.4348  1.2266  0.2335
-0.069   0.1286  0.6134
```
And then inverted with MATLAB:
```matlab
// https://www.mathworks.com/help/matlab/ref/inv.html
X = [1.3717 -0.649 -0.1154; -0.4348 1.2266 0.2335; -0.069 0.1286 0.6134]
Y = inv(X)
```
```
0.8757	0.4646	-0.0121
0.3038	1.0103	-0.3274
0.0348	-0.1596	1.6975
```



