float WorldSize = 10.0;
int NX = 64;
int NY = 64;
int ArraySize = NX * NY;
float DXY = WorldSize / NX;

float LX = WorldSize;
float LY = WorldSize;

int StateSize = 3;
float[][] State = new float[StateSize][ArraySize];
int StateVelU = 0;
int StateVelV = 1;
int StateSoot = 2;

int PixelsPerCell = 8;

int WindowWidth = PixelsPerCell * NX;
int WindowHeight = PixelsPerCell * NY;

PImage StateImage = createImage( NX, NY, RGB );

// Index an element of a grid in the state array
int IX( int i, int j ) {
    return ( i + NX*j );
}

float snoise( float x, float y ) {
   return (2.0 * noise( x, y )) - 1.0;
}

float Fractal(float x, float y) {
  float octave_gain = 0.5;
  float lacunarity = 1.67;
  float f = 0.0;
  float amp = 1.0;
  float px = x;
  float py = y;
  for (int oct = 0; oct < 4; ++oct) {
    f += amp * snoise(px, py);
    amp *= octave_gain;
    px *= lacunarity;
    py *= lacunarity;
  }
  return f;
}

void SetFieldToFractal(int field, float gain,
    float offset_x, float offset_y,
    float pattern_size) {
  noiseSeed( 0 );
  float spacing = DXY / pattern_size;
  for (int j = 0; j < NY; ++j) {
    float noise_pos_y = offset_y + spacing * float(j);    
    for (int i = 0; i < NX; ++i) {
      float noise_pos_x = offset_x + spacing * float(i);
      float d = gain * Fractal(noise_pos_x, noise_pos_y);
      State[field][IX(i,j)] = d;
    }
  }
}

// Draws a hatched grid of lines.
void SetFieldToHatch(int field, int line_width, int cell_width) {
  for (int j = 0; j < NY; ++j) {
    int j_line = (j + (line_width/2)) % cell_width;
    boolean j_in_line = j_line < line_width;
    for (int i = 0; i < NX; ++i) {
      int i_line = (i + (line_width/2)) % cell_width;
      boolean i_in_line = i_line < line_width;
      
      State[field][IX(i,j)] = i_in_line || j_in_line ? 1.0 : 0.0;
    }
  } 
}

void SetInitialState() {
  SetFieldToFractal(StateVelV, WorldSize * 0.25, 
      2341.17, 9911.44, WorldSize * 0.617);
  SetFieldToFractal(StateVelU, WorldSize * 0.25, 
      -8181.13, -1881.81, WorldSize * 0.617);
  SetFieldToHatch(StateSoot, 2, NX/8);
}

void setup() {
    SetInitialState();
    size(WindowWidth, WindowHeight);
    colorMode(RGB, 1.0);
    strokeWeight(0.5);
    textSize(24);
    noLoop();
}

// Display Gamma
float DisplayGamma = 2.2;

// Draw vector field into the image. The from_low and from_high 
// represent the numerical range that will be mapped into the 0-1
// space of red and green. A display gamma will then be applied.
void DrawVectorField(int field_r, int field_g, int field_b,
    float from_low_r, float from_low_g, float from_low_b,
    float from_high_r, float from_high_g, float from_high_b) {
  float r, g, b;
  StateImage.loadPixels();
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      r = field_r < 0 ? 0 : State[field_r][IX(i,j)];
      g = field_g < 0 ? 0 : State[field_g][IX(i,j)];
      b = field_b < 0 ? 0 : State[field_b][IX(i,j)];
      
      // remap.
      r = (r - from_low_r) / (from_high_r - from_low_r);
      g = (g - from_low_g) / (from_high_g - from_low_g);
      b = (b - from_low_b) / (from_high_b - from_low_b);
      
      // constrain.
      r = constrain(r, 0.0, 1.0);
      g = constrain(g, 0.0, 1.0);
      b = constrain(b, 0.0, 1.0);
      
      // display gamma.
      r = pow(r, DisplayGamma);
      g = pow(g, DisplayGamma);
      b = pow(b, DisplayGamma);
      
      // Set
      StateImage.pixels[IX(i,j)] = color(r, g, b);
    }
  }
  StateImage.updatePixels();
  image(StateImage, 0, 0, width, height);
}

void draw() {
  // Draw Velocity 
  DrawVectorField(StateVelU, StateVelV, StateSoot,
    -1.5, -1.5, 0.0,
    1.5, 1.5, 1.0);

  // Label.
  noStroke();
  fill(1.0, 1.0, 1.0, 0.75);
  rect(0.0, 0.0, width, 40.0);
  fill( 0.0 );
  text("Fluid Velocity & Soot : RG, B", 10, 30);
}

