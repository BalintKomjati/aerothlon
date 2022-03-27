
# Animated visualization of digital elevation data in R with satellite overlay and GPS tracklogs

<https://user-images.githubusercontent.com/47415815/139556600-15f502df-a3ab-4b67-8073-825972a3e752.mp4>

## Intro

[Aerothlon](https://www.aerothlon.com) is an alpine triathlon race of
mountain running, mountain biking and paragliding. I had the pleasure to
take part in the 2021 edition in Kleinarl (AT). This work visualizes the
efforts of all athletes at 350x speed.

## Give it a try

1.  Clone the repository into a new Rstudio project

2.  Install the following packages

``` r
  install.packages('tictoc') 
  install.packages('rayshader')
  install.packages('raster')
  install.packages('geoviz')
  install.packages('data.table')
```

3.  Install [ffmpeg](https://www.ffmpeg.org/) on your computer to render
    the video, and give R the ffmpeg’s install folder at the top of the
    the aerothlon.R script

``` r
ffmpeg_install_folder = "C:\\ffmpeg\\bin;"
```

4.  And then just source the script

``` r
source("aerothlon.R")
```

This will load the demo data from disk, generate the frames 1-by-1, and
render the video into the working directory as ‘animation.mp4’.

The script is loaded with comments so feel free to sneak in and adjust.

## Credits

This work utilizes the one-off functionalities of the awsome packages
[rayshader](https://www.rayshader.com/) &
[geoviz](https://github.com/neilcharles/geoviz)!

To large extend this visualization was inspired by (and built on top of)
the package
[rayshaderanimate](https://github.com/zappingseb/rayshaderanimate)
however i decided to write a dedicated animation function on my own due
to the massive changes needed in the function
rayshaderanimate::video\_animation\_rayshade() to fit my needs.

## Notes

You also need to register at <https://www.mapbox.com/> for a free api
key to download custom overlay images (a low quality demo overlay is
included)

The script above is designed to work with a list of tracklogs however
for this demo it is supplied with only one. Since it was a paragliding
competition the logs are in igc format. The script should work well with
gpx files also. I preformatted the igc\_demo file e.g. adjusted the
elevation not to cross the ground and reduced its resolution.
