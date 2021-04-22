# Cropper

A very quickly put-together (Linux-only) image cropping tool.

If you have a large number of images, and you want to crop out bits of those
images all with the same aspect ratio, this is the tool for you.

To compile, run `make release`.

There may be bugs; I won't fix them (sorry).

See the top of `main.c` for basic configuration.
Click to make a new crop. The first one will be stored as `[filename]-1.png`, then
`[filename]-2.png`, etc. You can press a number key to switch to that index.

Scroll to resize the cropping frame, and press r to reset it to the default.
Press space/n to move to the next image, and p to move to the previous one.

