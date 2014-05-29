float WaveSpeed = 0.5;

float WorldSize = 10.0;
int NX = 64;
int NY = 64;
int ArraySize = NX * NY;
float DXY = WorldSize / NX;

float LX = WorldSize;
float LY = WorldSize;
float LZ = WorldSize / 2.0;

int StateSize = 8;
float[][] State = new float[StateSize][ArraySize];
int StateHeight = 0;
int StateVel = 1;
int StateHeightPrev = 2;
int StateVelPrev = 3;
int StateVelStar = 4;
int StateAccelStar = 5;
int StateHeightStar = 6;
int StateJacobiTmp = 7;

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

void EnforceHeightBoundaryConditions( int io_h ) {
    EnforceNeumannBoundaryConditions(io_h);

    if ( InputActive ) {
        State[io_h][IX(InputIndexX, InputIndexY)] = InputHeight;
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

void SetInitialState() {
    noiseSeed( 0 );
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            float worldX = 2341.17 + DXY * ( float )i;
            float worldY = 9911.44 + DXY * ( float )j;

            float n = 0.5 * snoise( worldX * 0.0625, worldY * 0.0625 ) +
                0.4 * snoise( worldX * 0.125, worldY * 0.125 ) +
                0.3 * snoise( worldX * 0.25, worldY * 0.25 ) +
                0.2 * snoise( worldX * 0.5, worldY * 0.5 );
            n = 1.0 - abs(n);
            n = (2.0 * n) - 1.0;

            State[StateHeight][IX(i,j)] = n;

            State[StateVel][IX(i,j)] = 0.0;
        }
    }
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceNeumannBoundaryConditions( StateVel );
    StateCurrentTime = 0.0;

    CopyArray( StateHeight, StateHeightPrev );
    CopyArray( StateVel, StateVelPrev );
    InputActive = false;
}

void setup() {
    SetInitialState();

    size( WindowWidth, WindowHeight );

    colorMode( RGB, 1.0 );
    strokeWeight( 0.5 );
    textSize( 24 );
}

// Draw height field into the image.
void DrawHeightField( int i_field ) {
    float pixr, pixg, pixb;
    float d;
    color dark_blue = color(0.01, 0.01, 0.2);
    color light_blue = color(0.9, 0.9, 1.0);
    StateImage.loadPixels();
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            float d = State[i_field][IX(i,j)];

            // Renormalize
            d = constrain( (d + 1.0) / 2.0, 0.0, 1.0 );

            // Add contrast
            d = pow( d, 8 );

            // Interpolate color.

            StateImage.pixels[IX(i,j)] = lerpColor(dark_blue, light_blue, d);
        }
    }
    StateImage.updatePixels();
    image( StateImage, 0, 0, width, height );
}

//-*****************************************************************************
void draw() {
    background( 0.5 );

    DrawHeightField(StateHeight);

    // Label.
    fill( 1.0 );
    text("2D Wave Equation : Init Only", 10, 30);
}

// Reset function. If the key 'r' is released in the display,
// copy the initial state to the state.
void keyReleased() {
    if ( key == 114 ) {
        SetInitialState();
    }
}
