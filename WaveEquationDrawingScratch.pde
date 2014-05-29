
float schlick(float n2, float cos_theta) {
    float n1 = 1.0;
    float R0 = sq((n1 - n2)/(n1 + n2));
    return R0 + (1.0 - R0)*pow(1.0 - cos_theta, 5.0);
}

float kSpecular(PVector In, PVector Nn, PVector Ln, float m )
{
    PVector Vn = PVector.mult(In, -1.0);
    PVector H = PVector.add(Ln, Vn);
    H.normalize();
    float d = Nn.dot(H);
    d *= d;
    return pow( d, m/2 );
}

float kDiffuse(PVector Nn, PVector Ln)
{
    float d = Nn.dot(Ln);
    return constrain( d, 0, 1 );
}

//-*****************************************************************************
// DRAW HEIGHT FIELD
//-*****************************************************************************
void DrawHeightField( int i_field ) {
    float pixr, pixg, pixb;
    float d;
    StateImage.loadPixels();
    PVector dh_di = new PVector(0, 0, 0);
    PVector dh_dj = new PVector(0, 0, 0);
    PVector n = new PVector(0, 0, 0);
    PVector toSun = new PVector(1, 1, 2.5);
    toSun.normalize();
    PVector toMoon = new PVector(-1, -.5, 1.25);
    toMoon.normalize();
    PVector In = new PVector(0, 0, -1.0);
    float specular_intensity_r = 30.0;
    float specular_intensity_g = 30.0;
    float specular_intensity_b = 30.0;
    float diffuse_intensity_r = 0.05;
    float diffuse_intensity_g = 0.125;
    float diffuse_intensity_b = 0.0625;
    float sun_r = 1.0;
    float sun_g = 1.0;
    float sun_b = 1.0;
    float moon_r = 0.1;
    float moon_g = 0.1;
    float moon_b = 0.25;
    for (int j = 0; j < NY; ++j) {
        int j_neg = constrain(j-1, 0, NY-1);
        int j_pos = constrain(j+1, 0, NY-1);
        float dy;
        if (j == 0 || j == (NY-1)) {
            dy = DXY;
        } else {
            dy = 2.0 * DXY;
        }
        for (int i = 0; i < NX; ++i) {
            int i_neg = constrain(i-1, 0, NX-1);
            int i_pos = constrain(i+1, 0, NX-1);
            float dx;
            if (i == 0 || i == (NX-1)) {
                dx = DXY;
            } else {
                dx = 2.0 * DXY;
            }
            float h_left = State[i_field][IX(i_neg,j)];
            float h_right = State[i_field][IX(i_pos,j)];
            float h_down = State[i_field][IX(i,j_neg)];
            float h_up = State[i_field][IX(i,j_pos)];
            float h_cen = State[i_field][IX(i,j)];

            dh_di.set(dx, 0.0, h_right - h_left);
            dh_dj.set(0.0f, dy, h_up - h_down);
            n.set(dh_di.cross(dh_dj));
            n.normalize();

            float fresnel_d = schlick(1.31, n.z);
            float sun_spec = fresnel_d * kSpecular(In, n, toSun, 80);//sq(n.dot(toSun));
            float sun_spec_r = specular_intensity_r * sun_spec * sun_r;
            float sun_spec_g = specular_intensity_g * sun_spec * sun_g;
            float sun_spec_b = specular_intensity_b * sun_spec * sun_b;

            float sun_diff = kDiffuse(n, toSun);
            float sun_diff_r = diffuse_intensity_r * sun_diff * sun_r;
            float sun_diff_g = diffuse_intensity_g * sun_diff * sun_g;
            float sun_diff_b = diffuse_intensity_b * sun_diff * sun_b;


            float moon_spec = fresnel_d * kSpecular(In, n, toMoon, 50);//sq(n.dot(toSun));
            float moon_spec_r = specular_intensity_r * moon_spec * moon_r;
            float moon_spec_g = specular_intensity_g * moon_spec * moon_g;
            float moon_spec_b = specular_intensity_b * moon_spec * moon_b;

            float moon_diff = kDiffuse(n, toMoon);
            float moon_diff_r = diffuse_intensity_r * moon_diff * moon_r;
            float moon_diff_g = diffuse_intensity_g * moon_diff * moon_g;
            float moon_diff_b = diffuse_intensity_b * moon_diff * moon_b;

            //d = constrain(State[i_field][IX(i,j)], -1.0, 1.0);
            //d = 0.5 + (0.5 * d);

            //colorMode(HSB, 1);
            //float h = map(d, 0, 1, 0.675, .55);
            //float s = map(d, 0, 1, 0, 1);
            //float b = map(d, 0, 1, 1, .25);
            //pixr = 0.9 * ( 1.0 - d );
            //pixg = 0.9 * ( 1.0 - (d*d) );
            //pixb = 0.9 * ( 1.0 - (d*d*d) );
            pixr = (sun_spec_r + sun_diff_r + moon_spec_r + moon_diff_r);
            pixg = (sun_spec_g + sun_diff_g + moon_spec_g + moon_diff_g);
            pixb = (sun_spec_b + sun_diff_b + moon_spec_b + moon_diff_b);

            StateImage.pixels[IX(i,j)] = color( pixr, pixg, pixb );
        }
    }
    StateImage.updatePixels();
    image( StateImage, 0, 0, width, height );
}
