movie:
	ffmpeg -framerate 10 -i %05d.png -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 10 movie.mp4
	#mencoder mf://*.png -fps 10 -ovc x264 -o movie.mp4

