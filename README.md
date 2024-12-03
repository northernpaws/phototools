# Raw Image Developer

C++-based Raw image developer, built on algorithms and principals originally found in DCRaw and Darktable. It works incredibly fast and has slightly higher-quality than DCRaw due to my changes in making the pipelines float-based.

Takes a source RAW image from a camera (only tested with Fujifilm XTrans RAW files so far), and convert it into a color image in the correct colorspace.

## The Process

First, the program needs to take the 1-bit grayscale image produced by the camera and de-mosaic it.

![Raw image data](/assets/raw.png)

This image appears as grayscale because the actual color data is encoded as a grid of different color filters in a known pattern; in our case, this is the XTrans pattern.

> ![XTrans camera sensor pattern](https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/XTrans_matrix.png/220px-XTrans_matrix.png)
>
> https://en.wikipedia.org/wiki/Fujifilm_X-Trans_sensor 

De-mosaicing can be a tricky process, we need to make sure we're processing the source image in the correct orientation and alignment to the XTrans pattern. If we don't process it in the right orientation, or process the pattern with the wrong alignment, we get results like this:

![Misaligned XTrans sensor data](/assets/demosaiced_missaligned.png)
_(NOTE: Image downscaled to save bandwidth)_

This looks neat, but isn't what we need. What we want to do in this first pass is to assign all the RGB pixels to their respective color channels in the image, so when we zoom in we should see a pattern like the source XTrans sensor pattern.

If we look a little closer, we can see the alignment pattern in the image.

![Missaligned XTrans sensor pattern](/assets/missaligned_pattern.png)

As we can see, the pattern we're seeing doesn't look similar to the reference pattern above. After shifting the alignment of the demosaicing algorithm slightly, we get a much better result:

![Correctly aligned demosaiced sensor data](/assets/demosaiced_aligned.png)
_(NOTE: Image downscaled to save bandwidth)_

The image appears more muted than the source image, but this is because it comes out of the camera in a neutral colorspace, and we'll need to convert that to a display colorspace. We're also viewing it as a grid of color channels and not a fully composed image. 

If we look up close this time, we can see the correct pattern:

![Correctly aligned color data](/assets/aligned_pattern.png)

Now that we have the color data sorted into the correct channels, we can run the interpolation algorithm.

The interpolation algorithm looks at each pixel in the image, and depending on what color and position it is relative to the sensor array, uses surrounding colors to interpolation the two missing color components for the pixel.

On standard color sensors (Bayer) this is pretty straightforward, however, Fuji's X-Trans layout is significantly more complicated. We've adapted the algorithm from DCRaw to perform this, as the colorspace theory required to write it from scratch is definitely over my head.  

![Decoded RAW image](/assets/decoded.JPG)

Note that we've also applied some colorspace transformations into display colorspace; otherwise the image would be much more muted.

# Dependencies

> Re-configure the cmake project a second time after the first failure to find automatically downloaded dependencies.

## Rawspeed

  - OpenMP
    - macOS
      - `brew install llvm libomp`
      - use /usr/local/opt/llvm/bin/clang for toolchain compiler (or disable openmp parallelization support!)

# Credits

- [Floating Point Dcraw + DCB, AMaZE, LMMSE, AFD, & VCD; Part 4](https://ninedegreesbelow.com/photography/dcraw-float-source-license.html) by NineDegreesBelow
  - Gave me a reference for what a float-based pipeline for DCRaw _could_ look like.
- [Darktable](https://www.darktable.org)
  - Used to generate proofs of my source RAW files using different algorithms to compare my program output against.
  - Develops the RawSpeed RAW file decoder used by Darktable.
- [Ansel](https://ansel.photos/)
  - Referenced for some of it's workflow improvements over Darktable.
- [FastRawViewer](https://www.fastrawviewer.com)
  - Used to export reference images from my source RAWs.
  - Reliable reference since it also uses libraw and rawspeed under the hood.

# License 

The code in this repository is licensed under the [GNU General Public License v2.0](/LICENSE.md) due to licensing requirements from the adapted DCRaw code.
