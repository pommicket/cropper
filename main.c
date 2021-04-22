// directory to take png images from
#define INPUT_DIR  "i"
// directory to output cropped png images to
#define OUTPUT_DIR "o"
// aspect ratio = width / height
#define ASPECT_RATIO (2.0f/3.0f)
// output width in pixels. all output images be resized to this width.
#define OUTPUT_WIDTH 1550
// starting scale for frame
#define STARTING_SCALE 0.8f

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wconversion"
#endif
#include <SDL.h>
#include <dirent.h>
#define STBI_ONLY_PNG
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_STATIC
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_STATIC
#include "stb_image_resize.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_STATIC
#include "stb_image_write.h"
#if __GNUC__
#pragma GCC diagnostic pop
#endif
#include <stdbool.h>
#if _WIN32 || WINDOWS
#include <windows.h>
#endif

#include "base.h"
#include "util.c"
#include "gl.c"
#include "math.c"
#include "time.c"
#include "arr.c"

static void die(char const *fmt, ...) {
	char buf[256] = {0};
	
	va_list args;
	va_start(args, fmt);
	vsnprintf(buf, sizeof buf - 1, fmt, args);
	va_end(args);

	// show a message box, and if that fails, print it
	if (SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, "Error", buf, NULL) < 0) {
		debug_println("%s\n", buf);
	}

	exit(EXIT_FAILURE);
}

#if DEBUG
static void APIENTRY gl_message_callback(GLenum source, GLenum type, unsigned int id, GLenum severity, 
	GLsizei length, const char *message, const void *userParam) {
	(void)source; (void)type; (void)id; (void)length; (void)userParam;
	if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;
	debug_println("Message from OpenGL: %s.", message);
}
#endif


static void load_texture(char const *filename, int *width, int *height, unsigned char *image, GLuint tex) {
	int c;
	int w, h;
	unsigned char *data = stbi_load(filename, &w, &h, &c, 1);
	if (data) {
		memcpy(image, data, (size_t)w*(size_t)h);
		// flip vertically
		for (int y = 0; y < h/2; ++y) {
			int oy = h-1-y;
			for (int x = 0; x < w; ++x) {
				unsigned char *a = &data[y*w+x], *b = &data[oy*w+x];
				unsigned char tmp = *a;
				*a = *b;
				*b = tmp;
			}
		}
		glBindTexture(GL_TEXTURE_2D, tex);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_UNSIGNED_BYTE, data);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		*width = w;
		*height = h;
		stbi_image_free(data);
	} else {
		printf("Error loading %s\n", filename);
	}
}

static int qsort_strcmp(const void *a, const void *b) {
	return strcmp(*(const char *const *)a, *(const char *const *)b);
}

int main(int argc, char **argv) {
	mkdir(OUTPUT_DIR, 0755);
	DIR *input_dir = opendir(INPUT_DIR);
	static char const *filenames[10000];
	int nfiles = 0;
	{
		struct dirent *ent;
		while ((ent = readdir(input_dir))) {
			if (ent->d_type == DT_REG) {
				char const *name = ent->d_name;
				if (strchr(name, '.') && strcmp(strchr(name, '.'), ".png") == 0)
					filenames[nfiles++] = strdup(ent->d_name);
			}
		}
		qsort(filenames, (size_t)nfiles, sizeof *filenames, qsort_strcmp);
	}
	if (nfiles == 0) {
		printf("No files!\n");
		return -1;
	}
	
	int file_idx = 0;
	if (argc < 2) {
		FILE *fp = fopen("pos.txt", "r");
		if (fp) {
			fscanf(fp, "%d", &file_idx);
			fclose(fp);
		}
	} else {
		file_idx = atoi(argv[1]);
	}
	
	setbuf(stdout, NULL);
		
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window *window = SDL_CreateWindow("cropper", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
		1280, 720, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);
	
	int gl_version_major = 0, gl_version_minor = 0;
	struct {
		int major;
		int minor;
	} gl_versions[] = {
	// first try gl_versions[0], then gl_versions[1], etc.
		{4, 3}, // debug context support
		{3, 2}, // framebuffer support
		{3, 0}, // vao support
	};
	
	SDL_GLContext *glctx = NULL;
	
	for (size_t i = 0; i < arr_count(gl_versions); ++i) {
		gl_version_major = gl_versions[i].major;
		gl_version_minor = gl_versions[i].minor;
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, gl_version_major);
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, gl_version_minor);
	#if DEBUG
		SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, gl_version_major * 100 + gl_version_minor >= 403 ? SDL_GL_CONTEXT_DEBUG_FLAG : 0);
	#endif
		
		glctx = SDL_GL_CreateContext(window);
		if (glctx) {
		#if DEBUG
			print("Got GL %d.%d context.\n", gl_version_major, gl_version_minor);
		#endif
			break;
		} else {
			debug_println("Couldn't get GL %d.%d context.", gl_version_major, gl_version_minor);
		}
	}
	if (!glctx)
		die("%s\n(your GPU/graphics drivers might be too old)", SDL_GetError());
	gl_get_procs();
	
#if DEBUG
	if (gl_version_major * 100 + gl_version_minor >= 403) {
		GLint flags = 0;
		glGetIntegerv(GL_CONTEXT_FLAGS, &flags);
		glEnable(GL_DEBUG_OUTPUT);
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
		if (flags & GL_CONTEXT_FLAG_DEBUG_BIT) {
			// set up debug message callback
			glDebugMessageCallback(gl_message_callback, NULL);
			glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_TRUE);
			printf("Set up debug message callback\n");
		}
	}
#endif

	SDL_GL_SetSwapInterval(1);

	
	Uint32 last_frame = SDL_GetTicks();
	
	bool quit = false;
	
	GLuint shader = gl_create_program_from_files("v.glsl", "f.glsl");
	GLuint vbo = gl_gen_buffer(), vao = gl_gen_vertex_array();
	GLuint v_pos = gl_attrib_loc(shader, "v_pos");
	GLuint v_tex_coord = gl_attrib_loc(shader, "v_tex_coord");
	GLuint v_color = gl_attrib_loc(shader, "v_color");
	GLint u_transform = gl_uniform_loc(shader, "u_transform");
	GLint u_aspect_ratio = gl_uniform_loc(shader, "u_aspect_ratio");
	GLint u_use_texture = gl_uniform_loc(shader, "u_use_texture");
	GLint u_texture = gl_uniform_loc(shader, "u_texture");
	
	int tex_width = 0, tex_height = 0;
	unsigned char *tex_image = malloc(32ul<<20);
	GLuint texture = 0;
	glGenTextures(1, &texture);
	char path[1024];
	sprintf(path, "%s/%s", INPUT_DIR, filenames[file_idx]);
	load_texture(path, &tex_width, &tex_height, tex_image, texture);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	float frame_aspect_ratio = ASPECT_RATIO;
	float frame_scale = OUTPUT_WIDTH * STARTING_SCALE;
	
	unsigned char *output = malloc(32ul<<20);
	
	int out_idx = 1;
	float out_idx_animation = 0;
	int out_idx_animation_number = 0;
	
	while (!quit) {
		Uint32 this_frame = SDL_GetTicks();	
		float dt = (float)(this_frame - last_frame) * 0.001f;
		dt = minf(dt, 0.1f); // set a maximum frame time to avoid weird edge cases involving really long frames
		last_frame = this_frame;
		
		float window_width, window_height;
		{
			int w = 0, h = 0;
			SDL_GetWindowSize(window, &w, &h);
			window_width = (float)w;
			window_height = (float)h;
			glViewport(0, 0, w, h);
		}
		float aspect_ratio = window_width / window_height;
		
		bool clicked = false;
		
		
		Uint8 const *keys_down = SDL_GetKeyboardState(NULL); (void)keys_down;
		
		SDL_Event event;
		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			case SDL_QUIT:
				quit = true;
				break;
			case SDL_KEYDOWN:
				switch (event.key.keysym.sym) {
				case SDLK_r:
					frame_scale = STARTING_SCALE * OUTPUT_WIDTH;
					break;
				case SDLK_n:
				case SDLK_SPACE:
					++file_idx;
					file_idx %= nfiles;
					goto reload;
				case SDLK_p:
					file_idx = (file_idx + nfiles-1) % nfiles;
					goto reload;
				reload:
					sprintf(path, "%s/%s", INPUT_DIR, filenames[file_idx]);
					load_texture(path, &tex_width, &tex_height, tex_image, texture);
					out_idx = 1;
					break;
				case SDLK_1: out_idx = 1; break;
				case SDLK_2: out_idx = 2; break;
				case SDLK_3: out_idx = 3; break;
				case SDLK_4: out_idx = 4; break;
				case SDLK_5: out_idx = 5; break;
				case SDLK_6: out_idx = 6; break;
				case SDLK_7: out_idx = 7; break;
				case SDLK_8: out_idx = 8; break;
				case SDLK_9: out_idx = 9; break;
				}
				break;
			case SDL_MOUSEWHEEL: {
				float scale_speed = 30.0f;
				frame_scale += (float)event.wheel.y * scale_speed;
				if (frame_scale / frame_aspect_ratio > tex_height) frame_scale = tex_height * frame_aspect_ratio;
				if (frame_scale > tex_width) frame_scale = tex_width;
				if (frame_scale < 50) frame_scale = 50;
			} break;
			case SDL_MOUSEBUTTONDOWN: {
				clicked = true;
			} break;
			}
		}
		
		int imouse_x = 0, imouse_y = 0;
		SDL_GetMouseState(&imouse_x, &imouse_y);
		imouse_y = (int)window_height-1-imouse_y;
		float mouse_x = ((float)imouse_x / window_width * 2 - 1) * aspect_ratio;
		float mouse_y = (float)imouse_y / window_height * 2 - 1;
		
		
		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		
		typedef struct {
			v2 pos;
			v2 tex_coord;
			v4 color;
		} Vertex;
		Vertex *vertices = NULL;
		
		int iframe_w = (int)frame_scale;
		int iframe_h = (int)frame_scale / frame_aspect_ratio;
		int mouse_tex_x, mouse_tex_y;
		float frame_dx, frame_dy; // half-width and half-height of cropping frame
		
		{
			float x = (float)tex_width / (float)tex_height;
			float y = 1;
			
			if (x > aspect_ratio) {
				x = 1;
				y = (float)tex_height / (float)tex_width;
			}
			
			frame_dx = (float)iframe_w / (float)tex_width  * x;
			frame_dy = (float)iframe_h / (float)tex_height * y;
			
			if (mouse_x - frame_dx < -x) mouse_x = -x + frame_dx;
			if (mouse_x + frame_dx > +x) mouse_x = +x - frame_dx;
			if (mouse_y - frame_dy < -y) mouse_y = -y + frame_dy;
			if (mouse_y + frame_dy > +y) mouse_y = +y - frame_dy;
			
			mouse_tex_x = (int)((float)tex_width  * (mouse_x / x * 0.5f + 0.5f));
			mouse_tex_y = tex_height-1-(int)((float)tex_height * (mouse_y / y * 0.5f + 0.5f));
			int cx1 = mouse_tex_x - iframe_w / 2;
			//int cx2 = mouse_tex_x + iframe_w / 2;
			int cy1 = mouse_tex_y - iframe_h / 2;
			//int cy2 = mouse_tex_y + iframe_h / 2;
			if (cx1 < 0) cx1 = 0;
			if (cy1 < 0) cy1 = 0;
			if (cx1 + iframe_w > tex_width)  iframe_w = tex_width  - cx1;
			if (cy1 + iframe_h > tex_height) iframe_h = tex_height - cy1;
			if (clicked) {
				int out_w = OUTPUT_WIDTH;
				int out_h = (int)(out_w / ASPECT_RATIO);
				stbir_resize_uint8(&tex_image[cy1*tex_width+cx1], iframe_w, iframe_h, tex_width,
					output, out_w, out_h, out_w, 1);
				sprintf(path, "%s/%.*s-%d.png", OUTPUT_DIR, (int)strcspn(filenames[file_idx], "."), filenames[file_idx], out_idx);
				out_idx_animation = 1;
				out_idx_animation_number = out_idx;
				stbi_write_png(path, out_w, out_h, 1, output, out_w);
				++out_idx;
			}
		
			arr_add(vertices, ((Vertex){{-x, -y}, {0, 0}, {1,1,1,1}}));
			arr_add(vertices, ((Vertex){{+x, -y}, {1, 0}, {1,1,1,1}}));
			arr_add(vertices, ((Vertex){{-x, +y}, {0, 1}, {1,1,1,1}}));
			arr_add(vertices, ((Vertex){{+x, +y}, {1, 1}, {1,1,1,1}}));
		}
		
		GLushort *elements = NULL;
		arr_add(elements, 0);
		arr_add(elements, 1);
		arr_add(elements, 2);
		arr_add(elements, 1);
		arr_add(elements, 3);
		arr_add(elements, 2);
		
		glBindVertexArray(vao);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, arr_size_in_bytes(vertices), vertices, GL_STREAM_DRAW);
		arr_clear(vertices);
		glVertexAttribPointer(v_pos, 2, GL_FLOAT, 0, sizeof(Vertex), (void *)offsetof(Vertex, pos));
		glVertexAttribPointer(v_tex_coord, 2, GL_FLOAT, 0, sizeof(Vertex), (void *)offsetof(Vertex, tex_coord));
		glVertexAttribPointer(v_color, 4, GL_FLOAT, 0, sizeof(Vertex), (void *)offsetof(Vertex, color));
		glEnableVertexAttribArray(v_pos);
		glEnableVertexAttribArray(v_tex_coord);
		glEnableVertexAttribArray(v_color);
		
		glUseProgram(shader);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);
		glUniform1i(u_use_texture, 1);
		glUniform1i(u_texture, 0);
		m3 transform = m3_identity;
		m3_uniform(u_transform, &transform);
		glUniform1f(u_aspect_ratio, aspect_ratio);
		glBindVertexArray(vao);
		glDrawElements(GL_TRIANGLES, arr_len(elements), GL_UNSIGNED_SHORT, elements);
		arr_clear(elements);
		
		{
			v4 c = {0.5f,0,0.5f,0.2f};
			float x1 = mouse_x - frame_dx;
			float x2 = mouse_x + frame_dx;
			float y1 = mouse_y - frame_dy;
			float y2 = mouse_y + frame_dy;
			// show frame
			arr_add(vertices, ((Vertex){{x1, y1}, {0, 0}, c}));
			arr_add(vertices, ((Vertex){{x2, y1}, {0, 0}, c}));
			arr_add(vertices, ((Vertex){{x1, y2}, {0, 0}, c}));
			arr_add(vertices, ((Vertex){{x2, y2}, {0, 0}, c}));
			arr_add(elements, 0);
			arr_add(elements, 1);
			arr_add(elements, 2);
			arr_add(elements, 1);
			arr_add(elements, 3);
			arr_add(elements, 2);
			
			if (out_idx_animation > 0) {
				v4 c2 = {0, 1, 0, out_idx_animation};
				// show out_idx
				out_idx_animation -= 0.0016f;
				int n = out_idx_animation_number;
				for (int i = 0; i < n; ++i) {
					float dx = 0.03f;
					float dy = 0.03f;
					float x = -0.9f + dx * 2 * 1.2f * (float)i;
					float y = -0.9f;
					GLushort e = (GLushort)arr_len(vertices);
					arr_add(elements, e+0);
					arr_add(elements, e+1);
					arr_add(elements, e+2);
					arr_add(elements, e+1);
					arr_add(elements, e+3);
					arr_add(elements, e+2);
					
					arr_add(vertices, ((Vertex){{x-dx, y-dy}, {0, 0}, c2}));
					arr_add(vertices, ((Vertex){{x+dx, y-dy}, {0, 0}, c2}));
					arr_add(vertices, ((Vertex){{x-dx, y+dy}, {0, 0}, c2}));
					arr_add(vertices, ((Vertex){{x+dx, y+dy}, {0, 0}, c2}));
				}
			}
		}
		
		glBindVertexArray(vao);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, arr_size_in_bytes(vertices), vertices, GL_STREAM_DRAW);
		arr_clear(vertices);
		glVertexAttribPointer(v_pos, 2, GL_FLOAT, 0, sizeof(Vertex), (void *)offsetof(Vertex, pos));
		glVertexAttribPointer(v_color, 4, GL_FLOAT, 0, sizeof(Vertex), (void *)offsetof(Vertex, color));
		glEnableVertexAttribArray(v_pos);
		glEnableVertexAttribArray(v_color);
		
		glUseProgram(shader);
		glUniform1i(u_use_texture, 0);
		m3_uniform(u_transform, &transform);
		glUniform1f(u_aspect_ratio, aspect_ratio);
		glBindVertexArray(vao);
		glDrawElements(GL_TRIANGLES, arr_len(elements), GL_UNSIGNED_SHORT, elements);
		arr_clear(elements);
		
		
		SDL_GL_SwapWindow(window);
	}
	
	SDL_DestroyWindow(window);
	FILE *fp = fopen("pos.txt", "w");
	if (fp) {
		fprintf(fp, "%d", file_idx);
		fclose(fp);
	}
	return 0;
}
