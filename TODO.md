

* implement other demosaicing options (bayer and generic)
  * https://rawpedia.rawtherapee.com/Demosaicing
* lensfun for len corrections (used by afinity photo)
* LibExif for exif data (decouple from raw reading)
* LibOpenJpeg vs libjpeg-turbo
* https://www.littlecms.com/color-engine/ for color management (used by afinity photo)
  * cp "/System/Library/ColorSync/Profiles/sRGB Profile.icc" "sRGB Profile.icc"
* Fix vcpkg not using llvm so we can enable openmp support on the libraw package
* https://opencolorio.readthedocs.io/en/latest/quick_start/for_devs.html for color managemnet (exports in afinity photo as OCIO adjustment module)
