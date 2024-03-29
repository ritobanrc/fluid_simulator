// shader.vert
#version 450

layout(location=0) in vec3 a_position;
layout(location=1) in vec3 a_color;

layout(location=0) out vec3 v_color;

layout(set=0, binding=0) // 1.
uniform Uniforms {
    mat4 u_view_proj; // 2.
    float u_point_size;
};

void main() {
    v_color = a_color;
    gl_Position = u_view_proj * vec4(a_position, 1.0);
    gl_PointSize = u_point_size / gl_Position.w;
}
