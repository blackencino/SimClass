float WorldSize = 10.0;
int NX = 64;
int NY = 64;
int ArraySize = NX * NY;
float DXY = WorldSize / NX;

float LX = WorldSize;
float LY = WorldSize;

int StateSize = 2;
float[][] State = new float[StateSize][ArraySize];
int StateVelU = 0;
int StateVelV = 1;

float StateCurrentTime = 0.0;

int PixelsPerCell = 8;

int WindowWidth = PixelsPerCell * NX;
int WindowHeight = PixelsPerCell * NY;

boolean InputActive = false;
int InputIndexX = 0;
int InputIndexY = 0;
float InputHeight = 0;

PImage StateImage = createImage( NX, NY, RGB );

// Index an element of a grid in the state array
int IX( int i, int j ) {
    return ( i + NX*j );
}

float snoise( float x, float y ) {
   return ( 2.0 * noise( x, y )) - 1.0;
}

void EnforceDirichletBoundaryConditions( int io_a ) {
    for (int j = 0; j < NY; ++j) {
        if (j == 0 || j == (NY-1)) {
            for (int i = 0; i < NX; ++i) {
                State[io_a][IX(i,j)] = 0.0;
            }
        } else {
            State[io_a][IX(0,j)] = 0.0;
            State[io_a][IX(NX-1,j)] = 0.0;
        }
    }
}

void EnforceNeumannBoundaryConditions( int io_v ) {
    for (int j = 0; j < NY; ++j) {
        if (j == 0) {
            for (int i = 0; i < NX; ++i) {
                State[io_v][IX(i,0)] = State[io_v][IX(i,1)];
            }
        } else if (j == (NY-1)) {
            for (int i = 0; i < NX; ++i) {
                State[io_v][IX(i,NY-1)] = State[io_v][IX(i,NY-2)];
            }
        }

        State[io_v][IX(0,j)] = State[io_v][IX(1,j)];
        State[io_v][IX(NX-1,j)] = State[io_v][IX(NX-2,j)];
    }
}

void CopyArray( int i_src, int o_dst ) {
    for ( int i = 0; i < ArraySize; ++i ) {
        State[o_dst][i] = State[i_src][i];
    }
}

void FillArray( int o_a, float i_val ) {
    for ( int i = 0; i < ArraySize; ++i ) {
        State[o_a][i] = i_val;
    }
}

float Fractal(float x, float y) {
  float f = 0.5 * snoise( x * 0.0625, y * 0.0625 ) +
            0.4 * snoise( x * 0.125, y * 0.125 ) +
            0.3 * snoise( x * 0.25, y * 0.25 ) +
            0.2 * snoise( x * 0.5, y * 0.5 );
  return f;
}

void SetFieldToFractal(int field, float gain, 
    float offset_x, float offset_y,
    float spacing_x, float spacing_y) {
  noiseSeed( 0 );
  for (int j = 0; j < NY; ++j) {
    float noise_pos_y = offset_y + spacing_y * (float)j;    
    for (int i = 0; i < NX; ++i) {
      float noise_pos_x = offset_x + spacing_x * (float)i;
      State[field][IX(i,j)] = gain * Fractal(noise_pos_x, noise_pos_y);
    }
  }
}

void SetInitialState() {
  SetFieldToFractal(StateVelU, 1.0, 2341.17, 9911.44, DXY, DXY);
  SetFieldToFractal(StateVelV, 1.0, -8181.13, -1881.81, DXY, DXY);
}

void setup() {
    SetInitialState();
    size(WindowWidth, WindowHeight);
    colorMode(RGB, 1.0);
    strokeWeight(0.5);
    textSize(24);
}

// Display Gamma
float DisplayGamma = 2.2;

// Draw vector field into the image. The from_low and from_high 
// represent the numerical range that will be mapped into the 0-1
// space of red and green. A display gamma will then be applied.
void DrawVectorField(int field_u, int field_v, 
    float from_low_u, float from_high_u,
    float from_low_v, float from_high_v) {
  float u, v;
  StateImage.loadPixels();
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      u = State[field_u][IX(i,j)];
      v = State[field_v][IX(i,j)];
      
      // remap.
      u = (u - from_low_u) / (from_high_u - from_low_u);
      v = (v - from_low_v) / (from_high_v - from_low_v);
      
      // constrain.
      u = constrain(u, 0.0, 1.0);
      v = constrain(v, 0.0, 1.0);
      
      // display gamma.
      u = pow(u, DisplayGamma);
      v = pow(v, DisplayGamma);
      
      // Set
      StateImage.pixels[IX(i,j)] = color(u, v, 0.0);
    }
  }
  StateImage.updatePixels();
  image(StateImage, 0, 0, width, height);
}

// Draw vector field as line segments. We can optionally skip
// some of the points, and we can scale the length of the vectors
// accordingly.
void DrawVectorFieldAsLines(int field_u, int field_v,
    int step_u, int step_v,
    float gain_u, float gain_v) {
  float u, v;
  beginShape(LINES);
  for (int j = step_v / 2; j < NY; j += step_v) {
    float start_y = PixelsPerCell * (float)j;
    for (int i = step_u / 2; i < NX; i += step_u) {
      float start_x = PixelsPerCell * (float)i;
      
      u = gain_u * State[field_u][IX(i,j)];
      v = gain_v * State[field_v][IX(i,j)];
      
      float end_x = start_x + PixelsPerCell * u;
      float end_y = start_y + PixelsPerCell * v;
      vertex(start_x, start_y);
      vertex(end_x, end_y);
    }
  }
  endShape();
}

void draw() {
  background(0.5);
  // Draw Velocity as lines
  // Length gain is somewhat arbitrary.
  float gain = 20.0 * (1.0 / 24.0) / DXY;
  stroke(0.6, 0.2, 0.1);
  DrawVectorFieldAsLines(StateVelU, StateVelV, 2, 2, gain, gain);

  // Label.
  fill( 1.0 );
  text("Fluid Velocity : Velocity Lines", 10, 30);
}