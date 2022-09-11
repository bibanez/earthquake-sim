x11.x: main.cc
	g++ main.cc -o x11.out -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

web: main.cc resize_shell.html
	emcc -o earthquake.html main.cc -Os -Wall ~/Sources/raylib/src/libraylib.a -I ~/Sources/raylib/src/ -L ~/Sources/raylib/src/libraylib.a -s USE_GLFW=3 --shell-file resize_shell.html -DPLATFORM_WEB
	mv earthquake.html index.html

zip: index.html earthquake.wasm earthquake.js
	zip earthquake.zip index.html earthquake.wasm earthquake.js
