ffmpeg -r 20 -f image2 -s 1000x1000 -i frame_%d.png -vcodec libx264 -crf 15 -pix_fmt yuv420p output.mp4
