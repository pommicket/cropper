ADDITIONAL_FLAGS=$(CFLAGS)  -Wall -Wextra -Wshadow -Wpedantic -pedantic -std=gnu11 -Wno-unused-function
DEBUG_CFLAGS=$(ADDITIONAL_FLAGS) -O0 -g -DDEBUG=1 -fsanitize=undefined -fno-sanitize-recover=undefined
RELEASE_CFLAGS=$(ADDITIONAL_FLAGS) -O3 -s
LIBS=-lGL -lSDL2 -I/usr/include/SDL2 -lm
cropper: *.[ch]
	$(CC) $(DEBUG_CFLAGS) -o cropper main.c $(LIBS)
release: *.[ch]
	$(CC) $(RELEASE_CFLAGS) -o cropper main.c $(LIBS)
