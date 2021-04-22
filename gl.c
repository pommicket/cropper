#include "lib/glcorearb.h"

// macro trickery to avoid having to write everything twice
#define gl_for_each_proc(do)\
	do(DRAWARRAYS, DrawArrays)\
	do(DRAWELEMENTS, DrawElements)\
	do(GENTEXTURES, GenTextures)\
	do(DELETETEXTURES, DeleteTextures)\
	do(TEXIMAGE2D, TexImage2D)\
	do(BINDTEXTURE, BindTexture)\
	do(TEXPARAMETERI, TexParameteri)\
	do(GETERROR, GetError)\
	do(GETINTEGERV, GetIntegerv)\
	do(ENABLE, Enable)\
	do(DISABLE, Disable)\
	do(BLENDFUNC, BlendFunc)\
	do(DEPTHFUNC, DepthFunc)\
	do(VIEWPORT, Viewport)\
	do(CLEARCOLOR, ClearColor)\
	do(CLEAR, Clear)\
	do(FINISH, Finish)\
	do(PIXELSTOREI, PixelStorei)\
	do(CREATESHADER, CreateShader)\
	do(DELETESHADER, DeleteShader)\
	do(CREATEPROGRAM, CreateProgram)\
	do(SHADERSOURCE, ShaderSource)\
	do(GETSHADERIV, GetShaderiv)\
	do(GETSHADERINFOLOG, GetShaderInfoLog)\
	do(COMPILESHADER, CompileShader)\
	do(CREATEPROGRAM, CreateProgram)\
	do(DELETEPROGRAM, DeleteProgram)\
	do(ATTACHSHADER, AttachShader)\
	do(LINKPROGRAM, LinkProgram)\
	do(GETPROGRAMIV, GetProgramiv)\
	do(GETPROGRAMINFOLOG, GetProgramInfoLog)\
	do(USEPROGRAM, UseProgram)\
	do(GETATTRIBLOCATION, GetAttribLocation)\
	do(GETUNIFORMLOCATION, GetUniformLocation)\
	do(GENBUFFERS, GenBuffers)\
	do(DELETEBUFFERS, DeleteBuffers)\
	do(BINDBUFFER, BindBuffer)\
	do(BUFFERDATA, BufferData)\
	do(VERTEXATTRIBPOINTER, VertexAttribPointer)\
	do(VERTEXATTRIBIPOINTER, VertexAttribIPointer)\
	do(ENABLEVERTEXATTRIBARRAY, EnableVertexAttribArray)\
	do(DISABLEVERTEXATTRIBARRAY, DisableVertexAttribArray)\
	do(GENVERTEXARRAYS, GenVertexArrays)\
	do(DELETEVERTEXARRAYS, DeleteVertexArrays)\
	do(BINDVERTEXARRAY, BindVertexArray)\
	do(ACTIVETEXTURE, ActiveTexture)\
	do(UNIFORM1F, Uniform1f)\
	do(UNIFORM2F, Uniform2f)\
	do(UNIFORM3F, Uniform3f)\
	do(UNIFORM4F, Uniform4f)\
	do(UNIFORM1I, Uniform1i)\
	do(UNIFORM2I, Uniform2i)\
	do(UNIFORM3I, Uniform3i)\
	do(UNIFORM4I, Uniform4i)\
	do(UNIFORMMATRIX4FV, UniformMatrix2fv)\
	do(UNIFORMMATRIX4FV, UniformMatrix3fv)\
	do(UNIFORMMATRIX4FV, UniformMatrix4fv)\
	do(DEBUGMESSAGECALLBACK, DebugMessageCallback)\
	do(DEBUGMESSAGECONTROL, DebugMessageControl)\
	do(POLYGONMODE, PolygonMode)\
	do(GENFRAMEBUFFERS, GenFramebuffers)\
	do(BINDFRAMEBUFFER, BindFramebuffer)\
	do(FRAMEBUFFERTEXTURE2D, FramebufferTexture2D)\
	do(DELETEFRAMEBUFFERS, DeleteFramebuffers)\
	do(BLITFRAMEBUFFER, BlitFramebuffer)\
	do(DRAWBUFFER, DrawBuffer)\
	do(TEXIMAGE2DMULTISAMPLE, TexImage2DMultisample)\
	do(GENRENDERBUFFERS, GenRenderbuffers)\
	do(DELETERENDERBUFFERS, DeleteRenderbuffers)\
	do(BINDRENDERBUFFER, BindRenderbuffer)\
	do(RENDERBUFFERSTORAGE, RenderbufferStorage)\
	do(RENDERBUFFERSTORAGEMULTISAMPLE, RenderbufferStorageMultisample)\
	do(FRAMEBUFFERRENDERBUFFER, FramebufferRenderbuffer)\
	do(CHECKFRAMEBUFFERSTATUS, CheckFramebufferStatus)

#define gl_declare_proc(upper, lower) static PFNGL##upper##PROC gl##lower;
gl_for_each_proc(gl_declare_proc)
#undef gl_declare_proc

static void gl_get_procs(void) {
#if __GNUC__ && !__clang__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wpedantic"
#endif
	#define gl_get_proc(upper, lower) gl##lower = (PFNGL##upper##PROC)SDL_GL_GetProcAddress("gl" #lower);
	gl_for_each_proc(gl_get_proc)
	#undef gl_get_proc
#if __GNUC__ && !__clang__
	#pragma GCC diagnostic pop
#endif
}

#define BUFFER_BATCH_SIZE 8
#define VARRAY_BATCH_SIZE 8
static GLuint gl_buffers[BUFFER_BATCH_SIZE], gl_buffers_used = BUFFER_BATCH_SIZE;
static GLuint gl_varrays[VARRAY_BATCH_SIZE], gl_varrays_used = VARRAY_BATCH_SIZE;


static GLuint gl_gen_buffer(void) {
	if (gl_buffers_used >= BUFFER_BATCH_SIZE) {
		gl_buffers_used = 0;
		glGenBuffers(BUFFER_BATCH_SIZE, gl_buffers);
	}
	return gl_buffers[gl_buffers_used++];
}

static GLuint gl_gen_vertex_array(void) {
	if (gl_varrays_used >= VARRAY_BATCH_SIZE) {
		gl_varrays_used = 0;
		glGenVertexArrays(VARRAY_BATCH_SIZE, gl_varrays);
	}
	return gl_varrays[gl_varrays_used++];
}


// compile a GLSL shader
static GLuint gl_compile_shader(char const *code, GLenum shader_type) {
	GLuint shader = glCreateShader(shader_type);
	glShaderSource(shader, 1, &code, NULL);
	glCompileShader(shader);
	GLint status = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE) {
		char log[1024] = {0};
		glGetShaderInfoLog(shader, sizeof log - 1, NULL, log);
		debug_println("Error compiling shader: %s", log);
		return 0;
	}
	return shader;
}

// link together GL shaders
static GLuint gl_link_program(GLuint *shaders, size_t count) {
	GLuint program = glCreateProgram();
	if (program) {
		for (size_t i = 0; i < count; ++i) {
			if (!shaders[i]) {
				glDeleteProgram(program);
				return 0;
			}
			glAttachShader(program, shaders[i]);
		}
		glLinkProgram(program);
		GLint status = 0;
		glGetProgramiv(program, GL_LINK_STATUS, &status);
		if (status == GL_FALSE) {
			char log[1024] = {0};
			glGetProgramInfoLog(program, sizeof log - 1, NULL, log);
			debug_println("Error linking shaders: %s", log);
			glDeleteProgram(program);
			return 0;
		}
	}
	return program;
}

static GLuint gl_compile_and_link_shaders(char const *vshader_code, char const *fshader_code) {
	GLuint shaders[2];
	shaders[0] = gl_compile_shader(vshader_code, GL_VERTEX_SHADER);
	shaders[1] = gl_compile_shader(fshader_code, GL_FRAGMENT_SHADER);
	GLuint program = gl_link_program(shaders, 2);
	if (shaders[0]) glDeleteShader(shaders[0]);
	if (shaders[1]) glDeleteShader(shaders[1]);
	if (program) {
		debug_print("Successfully linked program %u.\n", program);
	}
	return program;
}

static GLuint gl_create_program_from_files(char const *vshader_file, char const *fshader_file) {
	char *vcode = read_file_contents(vshader_file);
	char *fcode = read_file_contents(fshader_file);
	GLuint program = 0;
	if (vcode && fcode) {
		program = gl_compile_and_link_shaders(vcode, fcode);
	} else {
		debug_print("WARNING: Couldn't read file: %s.\n", vcode ? fshader_file : vshader_file);
	}
	free(vcode); free(fcode);
	return program;
}

static GLuint gl_attrib_loc(GLuint program, char const *attrib) {
	GLint loc = glGetAttribLocation(program, attrib);
	if (loc == -1) {
		debug_print("Couldn't find vertex attribute %s.\n", attrib);
		return 0;
	}
	return (GLuint)loc;
}

static GLint gl_uniform_loc(GLuint program, char const *uniform) {
	GLint loc = glGetUniformLocation(program, uniform);
	if (loc == -1) {
		debug_print("Couldn't find uniform: %s.\n", uniform);
		return -1;
	}
	return loc;
}
