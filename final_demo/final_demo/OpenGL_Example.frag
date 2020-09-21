#version 330

in vec3 frag_color;
in vec3 frag_normal;
in vec4 position;
in vec2 uv;

uniform vec4 light_position;
uniform vec4 world_origin;
uniform vec3 light_color_1;
uniform vec3 light_color_2;
uniform sampler2D texture_Colors;
uniform bool sky;

out vec3 out_color;

void main(void) {
	vec3 color;
	vec3 frag_normals = normalize(frag_normal);
	vec4 U = normalize(position - world_origin);
	vec3 L = vec3(light_position - position);
	float lightDistance = length(L);
	L = L / lightDistance;

	if (!sky){
		color = texture2D(texture_Colors, uv ).rgb * frag_color;
		out_color = light_color_1 * color * (0.5 + dot(0.5, dot(frag_normals, vec3(U)))) + light_color_2 * color * max(0, dot(frag_normals, L));
	}
	else
		out_color = texture2D(texture_Colors, uv ).rgb * frag_color;
}