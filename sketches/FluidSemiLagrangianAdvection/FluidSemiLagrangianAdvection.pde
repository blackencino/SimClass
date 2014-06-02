float WorldSize = 10.0;
int NX = 64;
int NY = 64;
int ArraySize = NX * NY;
float DXY = WorldSize / NX;

float LX = WorldSize;
float LY = WorldSize;

int StateSize = 7;
float[][] State = new float[StateSize][ArraySize];
int StateVelU = 0;
int StatePrevVelU = 1;
int StateVelV = 2;
int StatePrevVelV = 3;
int StateSoot = 4;
int StatePrevSoot = 5;
int StateTmp = 6;

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

void SetVelocityToCurlOfQ(int U, int V, int Q) {
  for (int j = 0; j < NY; ++j) {
    int j_lo = constrain(j-1, 0, NY-1);
    int j_hi = constrain(j+1, 0, NY-1);
    float dy = DXY * float(j_hi - j_lo);
    for (int i = 0; i < NX; ++i) {
      int i_lo = constrain(i-1, 0, NX-1);
      int i_hi = constrain(i+1, 0, NX-1);
      float dx = DXY * float(i_hi - i_lo);
      
      float dQdx = (State[Q][IX(i_hi,j)] - State[Q][IX(i_lo,j)])/dx;
      float dQdy = (State[Q][IX(i,j_hi)] - State[Q][IX(i,j_lo)])/dy;
      
      State[U][IX(i,j)] = dQdy;
      State[V][IX(i,j)] = -dQdx;
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

void CopyField(int src, int dst) {
  for (int i = 0; i < ArraySize; ++i) {
    State[dst][i] = State[src][i];
  }  
}

void FillField(int field, float value) {
  for (int i = 0; i < ArraySize; ++i) {
    State[field][i] = value;
  }
}

void SetInitialState() {
  SetFieldToFractal(StateTmp, 5.0, 2341.17, 9911.44, WorldSize * 1.7);
  SetVelocityToCurlOfQ(StateVelU, StateVelV, StateTmp);
  //FillField(StateVelU, 1.0);
  //FillField(StateVelV, 0.0);
  //SetFieldToFractal(StateVelV, WorldSize * 0.25, 
  //    2341.17, 9911.44, WorldSize * 0.617);
  //SetFieldToFractal(StateVelU, WorldSize * 0.25, 
  //    -8181.13, -1881.81, WorldSize * 0.617);
  SetFieldToHatch(StateSoot, 2, NX/8);
  CopyField(StateVelU, StatePrevVelU);
  CopyField(StateVelV, StatePrevVelV);
  CopyField(StateSoot, StatePrevSoot);
}

void setup() {
    SetInitialState();
    size(WindowWidth, WindowHeight);
    colorMode(RGB, 1.0);
    strokeWeight(1.0);
    textSize(24);
}

void SwapSoot() {
  int tmp = StatePrevSoot;
  StatePrevSoot = StateSoot;
  StateSoot = tmp;
}

void SwapVel() {
  int tmp = StatePrevVelU;
  StatePrevVelU = StateVelU;
  StateVelU = tmp;
  
  tmp = StatePrevVelV;
  StatePrevVelV = StateVelV;
  StateVelV = tmp;
}

// Bi-Linear Resampling of a field Q at some point fi, fj
float BiLinearResample(int Q, float fi, float fj) {
  int i_lo = floor(fi);
  float s = fi - float(i_lo);
  i_lo = constrain(i_lo, 0, NX-1);
  int i_hi = min(i_lo + 1, NX-1);
  
  int j_lo = floor(fj);
  float t = fj - float(j_lo);
  j_lo = constrain(j_lo, 0, NY-1);
  int j_hi = min(j_lo + 1, NY-1);
  
  float q00 = State[Q][IX(i_lo,j_lo)];
  float q10 = State[Q][IX(i_hi,j_lo)];
  float q01 = State[Q][IX(i_lo,j_hi)];
  float q11 = State[Q][IX(i_hi,j_hi)];
  return lerp(lerp(q00, q10, s), lerp(q01, q11, s), t);
}

// Semi-Lagrangian Advection of Q by U & V
void SemiLagrangianAdvect(int new_Q, int old_Q, int U, int V, float dt) {
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      float vel_u = State[U][IX(i,j)];
      float vel_v = State[V][IX(i,j)];
      
      float cell_vel_u = vel_u / DXY;
      float cell_vel_v = vel_v / DXY;
      
      float fi = float(i) - (dt * cell_vel_u);
      float fj = float(j) - (dt * cell_vel_v);
      
      State[new_Q][IX(i,j)] = BiLinearResample(old_Q, fi, fj);
    }  
  } 
}

void TimeStep(float dt) {
  SwapSoot();
  SwapVel();
  SemiLagrangianAdvect(StateSoot, StatePrevSoot, 
      StatePrevVelU, StatePrevVelV, dt);
  SemiLagrangianAdvect(StateVelU, StatePrevVelU, 
      StatePrevVelU, StatePrevVelV, dt);
  SemiLagrangianAdvect(StateVelV, StatePrevVelV, 
      StatePrevVelU, StatePrevVelV, dt);
}

// Display Gamma
float DisplayGamma = 2.2;

// Draw velocity as line segments. We can optionally skip
// some of the points, and we can scale the length of the vectors
// accordingly.
void DrawVelocityAsLines(int U, int V,
    int step_u, int step_v,
    float dt) {
  beginShape(LINES);
  for (int j = step_v / 2; j < NY; j += step_v) {
    float pixel_start_y = PixelsPerCell * float(j);
    for (int i = step_u / 2; i < NX; i += step_u) {
      float pixel_start_x = PixelsPerCell * float(i);
      
      float vel_u = State[U][IX(i,j)];
      float vel_v = State[V][IX(i,j)];
      float pixel_vel_u = PixelsPerCell * State[U][IX(i,j)] / DXY;
      float pixel_vel_v = PixelsPerCell * State[V][IX(i,j)] / DXY;
      
      float pixel_end_x = pixel_start_x + (dt * pixel_vel_u);
      float pixel_end_y = pixel_start_y + (dt * pixel_vel_v);
      vertex(pixel_start_x, pixel_start_y);
      vertex(pixel_end_x, pixel_end_y);
    }
  }
  endShape();
}

void DrawScalarField(int field, 
    float from_low, float from_high, float exponent,
    color color_low, color color_high) {
  float d;
  StateImage.loadPixels();
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      d = State[field][IX(i,j)];
      
      d = (d - from_low)/(from_high - from_low);
      d = constrain(d, 0.0, 1.0);
      d = pow(d, exponent);
      
      StateImage.pixels[IX(i,j)] = lerpColor(color_low, color_high, d);
    }
  }  
  StateImage.updatePixels();
  image(StateImage, 0, 0, width, height);
}

void draw() {
  TimeStep(1.0/24.0);
  
  background(0.5);
  // Draw soot.
  color soot_high = color(0.3, 0.0, 0.3);
  color soot_low = color(1.0, 1.0, 1.0);
  DrawScalarField(StateSoot, 0.0, 1.0, 3.0, soot_low, soot_high);
  
  // Draw Velocity as lines
  // Length gain is somewhat arbitrary.
  float dt = (1.0 / 24.0);
  stroke(1.0, 0.1, 0.1);
  DrawVelocityAsLines(StateVelU, StateVelV, 2, 2, 10.0 * dt);

  // Label.
  noStroke();
  fill(1.0, 1.0, 1.0, 0.75);
  rect(0.0, 0.0, width, 40.0);
  fill( 0.0 );
  text("Semi-Lagrangian Advection", 10, 30);
}

// Reset function. If the key 'r' is released in the display,
// copy the initial state to the state.
void keyReleased() {
    if ( key == 114 ) {
        SetInitialState();
    }
}
