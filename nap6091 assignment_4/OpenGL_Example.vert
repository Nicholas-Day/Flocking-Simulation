#version 330

uniform mat4 M;
uniform mat4 MVP;
uniform mat3 normalMatrix;

in  vec3 in_position;
in  vec3 in_color;
in  vec4 in_normals;
in  vec2 in_uv;

out vec3 frag_color;
out vec3 frag_normal;
out vec4 position;
out vec2 uv;

void main(void) {

    position = M * vec4(in_position.x, in_position.y, in_position.z, 1.0);
	gl_Position = MVP * vec4(in_position.x, in_position.y, in_position.z, 1.0);

	//mat3 normMatrix = transpose(inverse(mat3(normalMatrix)));

    frag_color = in_color;
	frag_normal = normalize(normalMatrix * vec3(in_normals));
	uv = in_uv;
}