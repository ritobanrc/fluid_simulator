#version 450

layout(location=0) in vec3 v_color;
layout(location=0) out vec4 f_color;

void main() {
    vec2 point = gl_PointCoord - vec2(0.5, 0.5);
    float radius = length(point);
    float alpha = 2. * clamp(0.5 - radius, 0., 1.);
    f_color = vec4(v_color, alpha);
}
