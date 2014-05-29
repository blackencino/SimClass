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

bool DisplayInstructions = true;

PImage StateImage = createImage( NX, NY, RGB );

// 0 for first-order,
// 1 for RK2
// 2 for RK4
int TimeStepMethod = 2;

// Index an element of a grid in the state array
int IX( int i, int j ) {
    return ( i + NX*j );
}

float snoise( float x, float y ) {
   return ( 2.0 * noise( x, y )) - 1.0;
}

float anoise( float x, float y ) {
   return ( -2.0 * abs( snoise( x, y )) ) + 1.0;
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
    DisplayInstructions = true;
}

void setup() {
    SetInitialState();

    size( WindowWidth, WindowHeight );

    colorMode( RGB, 1.0 );
    strokeWeight( 0.5 );
    textSize( 24 );
}

void SwapHeight() {
    int tmp = StateHeight;
    StateHeight = StateHeightPrev;
    StateHeightPrev = tmp;
}

void SwapVel() {
    int tmp = StateVel;
    StateVel = StateVelPrev;
    StateVelPrev = tmp;
}

void SwapState() {
    SwapHeight();
    SwapVel();
}

void GetInput() {
    InputActive = false;
    if ( mousePressed && mouseButton == LEFT ) {
        float mouseCellX = mouseX / ( ( float )PixelsPerCell );
        float mouseCellY = mouseY / ( ( float )PixelsPerCell );

        int iX = ( int )floor( mouseCellX + 0.5 );
        int iY = ( int )floor( mouseCellY + 0.5 );
        if ( iX > 0 && iX < NX-1 && iY > 0 && iY < NY-1 ) {
            InputIndexX = iX;
            InputIndexY = iY;
            InputHeight = 1.5;
            InputActive = true;
            DisplayInstructions = false;
        }
    }
}

// Estimate height star
void EstimateHeightStar( float i_dt ) {
    for ( int i = 0; i < ArraySize; ++i ) {
        State[StateHeightStar][i] = State[StateHeightPrev][i] +
                ( i_dt * State[StateVelStar][i] );
    }
    EnforceHeightBoundaryConditions( StateHeightStar );
}

// Estimate vel star
void EstimateVelStar( float i_dt ) {
    for ( int i = 0; i < ArraySize; ++i ) {
        State[StateVelStar][i] = State[StateVelPrev][i] +
                ( i_dt * State[StateAccelStar][i] );
    }
    EnforceNeumannBoundaryConditions( StateVelStar );
}

// Jacobi iteration to get temp acceleration
void JacobiIterationAccel( int i_aOld, int o_aNew, int i_hStar, float i_dt ) {
    float kappa = sq( WaveSpeed ) * sq( i_dt ) / sq( DXY );
    float gamma = sq( WaveSpeed ) / sq( DXY );

    for (int j = 1; j < NY-1; ++j) {
        for (int i = 1; i < NX-1; ++i) {
            float a_left = State[i_aOld][IX(i-1,j)];
            float a_right = State[i_aOld][IX(i+1,j)];
            float a_down = State[i_aOld][IX(i,j-1)];
            float a_up = State[i_aOld][IX(i,j+1)];

            float h_star_left = State[i_hStar][IX(i-1,j)];
            float h_star_right = State[i_hStar][IX(i+1,j)];
            float h_star_down = State[i_hStar][IX(i,j-1)];
            float h_star_up = State[i_hStar][IX(i,j+1)];
            float h_star_cen = State[i_hStar][IX(i,j)];

            float b = gamma *
                (h_star_left + h_star_right + h_star_down + h_star_up -
                 (4.0 * h_star_cen));

            float c = kappa * (a_left + a_right + a_down + a_up);

            State[o_aNew][IX(i,j)] = (b + c) / (1.0 + kappa);
        }
    }

    EnforceDirichletBoundaryConditions( o_aNew );
}

// Solve for acceleration.
void JacobiSolveAccel( int i_hStar, float i_dt ) {
    // Initialize acceleration to zero.
    FillArray( StateAccelStar, 0.0 );

    // Solve from StateJacobiTmp into StateAccel
    for ( int iter = 0; iter < 20; ++iter ) {
        int tmp = StateAccelStar;
        StateAccelStar = StateJacobiTmp;
        StateJacobiTmp = tmp;

        JacobiIterationAccel( StateJacobiTmp, StateAccelStar, i_hStar, i_dt );
    }
}

void EstimateAccelStar( float i_dt ) {
    JacobiSolveAccel( StateHeightStar, i_dt );
}

// Accumulate estimate
void AccumulateEstimate( float i_dt ) {
    for ( int i = 0; i < ArraySize; ++i ) {
        State[StateHeight][i] += i_dt * State[StateVelStar][i];
        State[StateVel][i] += i_dt * State[StateAccelStar][i];
    }
}

// First-Order Symplectic Time Step function.
void TimeStepFirstOrder( float i_dt ) {
    // Swap state
    SwapState();

    // Initialize estimate. This just amounts to copying
    // The previous values into the current values.
    CopyArray( StateHeightPrev, StateHeight );
    CopyArray( StateVelPrev, StateVel );

    // 1
    CopyArray( StateVel, StateVelStar );
    EstimateHeightStar( i_dt );
    EstimateAccelStar( i_dt );
    AccumulateEstimate( i_dt );

    // Final boundary conditions on height and vel
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceNeumannBoundaryConditions( StateVel );

    // Update current time.
    StateCurrentTime += i_dt;
}

// Second-Order Runge-Kutta Time Step function.
void TimeStepRK2( float i_dt ) {
    // Swap state
    SwapState();

    // Initialize estimate. This just amounts to copying
    // The previous values into the current values.
    CopyArray( StateHeightPrev, StateHeight );
    CopyArray( StateVelPrev, StateVel );

    // 1
    CopyArray( StateVel, StateVelStar );
    EstimateHeightStar( i_dt );
    EstimateAccelStar( i_dt );
    AccumulateEstimate( i_dt / 2.0 );

    // 2
    EstimateVelStar( i_dt );
    EstimateHeightStar( i_dt );
    EstimateAccelStar( i_dt );
    AccumulateEstimate( i_dt / 2.0 );

    // Final boundary conditions on height and vel
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceNeumannBoundaryConditions( StateVel );

    // Update current time.
    StateCurrentTime += i_dt;
}

// Fourth-Order Runge-Kutta Time Step function.
void TimeStepRK4( float i_dt ) {
    // Swap state
    SwapState();

    // Initialize estimate. This just amounts to copying
    // The previous values into the current values.
    CopyArray( StateHeightPrev, StateHeight );
    CopyArray( StateVelPrev, StateVel );

    // 1
    CopyArray( StateVel, StateVelStar );
    EstimateHeightStar( i_dt );
    EstimateAccelStar( i_dt );
    AccumulateEstimate( i_dt / 6.0 );

    // 2
    EstimateVelStar( i_dt / 2.0 );
    EstimateHeightStar( i_dt / 2.0 );
    EstimateAccelStar( i_dt / 2.0 );
    AccumulateEstimate( i_dt / 3.0 );

    // 3
    EstimateVelStar( i_dt / 2.0 );
    EstimateHeightStar( i_dt / 2.0 );
    EstimateAccelStar( i_dt / 2.0 );
    AccumulateEstimate( i_dt / 3.0 );

    // 4
    EstimateVelStar( i_dt );
    EstimateHeightStar( i_dt );
    EstimateAccelStar( i_dt );
    AccumulateEstimate( i_dt / 6.0 );

    // Final boundary conditions on height and vel
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceNeumannBoundaryConditions( StateVel );

    // Update current time.
    StateCurrentTime += i_dt;
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

    GetInput();
    float dt = 1.0 / 24.0;
    if (TimeStepMethod == 0) {
        TimeStepFirstOrder(dt);
    } else if (TimeStepMethod == 1) {
        TimeStepRK2(dt);
    } else {
        TimeStepRK4(dt);
    }
    DrawHeightField(StateHeight);

    // Label.
    fill( 1.0 );
    if (TimeStepMethod == 0) {
        text("2D Wave Equation : Jacobi, First Order", 10, 30);
    } else if (TimeStepMethod == 1) {
        text("2D Wave Equation : Jacobi, RK2", 10, 30);
    } else {
        text("2D Wave Equation : Jacobi, RK4", 10, 30);
    }

    // Instructions
    if (DisplayInstructions) {
        text("Draw in me!", 200, 270);
    }
}

// Reset function. If the key 'r' is released in the display,
// copy the initial state to the state.
void keyReleased() {
    if ( key == 114 ) {
        SetInitialState();
    }
    if ( key == 116 ) {
        TimeStepMethod = (TimeStepMethod+1)%3;
    }
}
