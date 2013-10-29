float WaveSpeed = 0.5;

float WorldSize = 10.0;
int ArraySize = 128;
float DX = WorldSize / ArraySize;

float LX = WorldSize;
float LY = WorldSize / 2.0;

int StateSize = 9;
float[][] State = new float[StateSize][ArraySize];
int StateHeight = 0;
int StateVel = 1;
int StateHeightPrev = 2;
int StateVelPrev = 3;
int StateVelStar = 4;
int StateAccelStar = 5;
int StateHeightTmp = 6;
int StateJacobiTmp1 = 7;
int StateJacobiTmp2 = 8;

float StateCurrentTime = 0.0;

int PixelsPerCell = 4;

int WindowWidth = PixelsPerCell * ArraySize;
int WindowHeight = WindowWidth / 2;

float snoise( float x )
{
   return ( 2.0 * noise( x )) - 1.0;
}

float anoise( float x )
{
   return ( -2.0 * abs( snoise( x )) ) + 1.0;
}

void EnforceBoundaryConditions( int io_a )
{
    State[io_a][0] = ( 2.0 * State[io_a][1] ) - State[io_a][2];
    State[io_a][ArraySize-1] = ( 2.0 * State[io_a][ArraySize-2] ) - State[io_a][ArraySize-3];
}

void SetInitialState()
{
    noiseSeed( 0 );
    for ( int i = 0; i < ArraySize; ++i )
    {
        float worldX = 2341.17 + DX * ( float )i;
        State[StateHeight][i] = 0.5 * anoise( worldX * 0.0625 ) +
                         0.4 * anoise( worldX * 0.125 ) +
                         0.3 * anoise( worldX * 0.25 ) +
                         0.2 * anoise( worldX * 0.5 );
        State[StateVel][i] = 0.0;
    }
    EnforceBoundaryConditions( StateHeight );
    EnforceBoundaryConditions( StateVel );
    StateCurrentTime = 0.0;
}

void setup()
{
    SetInitialState();
    
    size( WindowWidth, WindowHeight );
    
    colorMode( RGB, 1.0 );
    strokeWeight( 0.5 );
    textSize( 24 );
}

void SwapHeight()
{
    int tmp = StateHeight;
    StateHeight = StateHeightPrev;
    StateHeightPrev = tmp;
}

void SwapVel()
{
    int tmp = StateVel;
    StateVel = StateVelPrev;
    StateVelPrev = tmp;
}

void SwapState()
{
    SwapHeight();
    SwapVel();
}

void CopyArray( int i_src, int o_dst )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[o_dst][i] = State[i_src][i];
    }
}

void GetInput()
{
    if ( mousePressed && mouseButton == LEFT )
    {
        float mouseCellX = mouseX / ( ( float )PixelsPerCell );
        float mouseCellY = ( (height/2) - mouseY ) / ( ( float )PixelsPerCell );
        float simY = mouseCellY * DX;
        
        int iX = ( int )floor( mouseCellX + 0.5 );
        if ( iX > 0 && iX < ArraySize-1 )
        {
            State[StateHeight][iX] = simY;
            State[StateVel][iX] = 0.0;
        }
    }
}

// Estimate temp height
void EstimateTempHeight( float i_dt )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[StateHeightTmp][i] = State[StateHeightPrev][i] + 
                ( i_dt * State[StateVelStar][i] );
    }
    EnforceBoundaryConditions( StateHeightTmp );
}

// Estimate vel star
void EstimateVelStar( float i_dt )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[StateVelStar][i] = State[StateVelPrev][i] + 
                ( i_dt * State[StateAccelStar][i] );
    }
    EnforceBoundaryConditions( StateVelStar );
}

// Accumulate estimate
void AccumulateEstimate( float i_dt )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[StateHeight][i] += i_dt * State[StateVelStar][i];
        State[StateVel][i] += i_dt * State[StateAccelStar][i];
    }
}

// Jacobi iteration to get temp height.
void JacobiIterationHeight( int i_h, int i_tmp, int i_v, float i_dt )
{
    float s = sq( WaveSpeed ) * sq( i_dt ) / sq( DX );
    
    int tmp = i_h;
    i_h = i_tmp;
    i_tmp = tmp;
    
    for ( int i = 1; i < ArraySize-1; ++i )
    {
        float hLeft = State[i_tmp][i-1];
        float hCen = State[i_tmp][i];
        float hRight = State[i_tmp][i+1];
        float vCen = State[i_v][i];
        
        float numer = hCen + ( vCen * i_dt ) + ( s * ( hLeft + hRight ) );
        float denom = 1.0 + ( 2.0 * s );
        
        State[i_h][i] = numer / denom;
    }
    
    EnforceBoundaryConditions( i_h );
}

void JacobiSolveHeight( int i_h, int i_tmp, int i_v, float i_dt )
{
    for ( int iter = 0; iter < 10; ++iter )
    {
        JacobiIterationHeight( i_h, i_tmp, i_v, i_dt );
    }
}

// Acceleration from height. This uses the spatial second derivative.
// Note that we're iterating only over the central points, not the edges,
// which are handled by the boundary condition.
void A_from_H( int i_h, int i_v, float i_dt )
{
    // Copy h to h_tmp
    CopyArray( i_h, StateJacobiTmp1 );
    
    // Jacobi-solve for new height that satisfies implicit equation.
    JacobiSolveHeight( StateJacobiTmp1, StateJacobiTmp2, i_v, i_dt );
  
    for ( int i = 1; i < ArraySize-1; ++i )
    {
        float hLeft = State[StateJacobiTmp1][i-1];
        float hCen = State[StateJacobiTmp1][i];
        float hRight = State[StateJacobiTmp1][i+1];

        float d2h_dx2 = ( hRight + hLeft - ( 2.0*hCen ) ) / sq( DX );

        State[StateAccelStar][i] = sq( WaveSpeed ) * d2h_dx2;
    }
    
    EnforceBoundaryConditions( StateAccelStar );
}

// Time Step function.
void TimeStep( float i_dt )
{
    // Swap state
    SwapState();
    
    // Initialize estimate. This just amounts to copying
    // The previous values into the current values.
    CopyArray( StateHeightPrev, StateHeight );
    CopyArray( StateVelPrev, StateVel );
    
    // Vstar1, Astar1
    CopyArray( StateVel, StateVelStar );
    A_from_H( StateHeight, StateVelStar, i_dt );
    // Accumulate
    AccumulateEstimate( i_dt / 6.0 );
    
    // Height Temp 2
    EstimateTempHeight( i_dt / 2.0 );
    // Vstar2, Astar2
    EstimateVelStar( i_dt / 2.0 );
    A_from_H( StateHeightTmp, StateVelStar, i_dt / 2.0 );
    // Accumulate
    AccumulateEstimate( i_dt / 3.0 );
    
    // Height Temp 3
    EstimateTempHeight( i_dt / 2.0 );
    // Vstar3, Astar3
    EstimateVelStar( i_dt / 2.0 );
    A_from_H( StateHeightTmp, StateVelStar, i_dt / 2.0 );
    // Accumulate
    AccumulateEstimate( i_dt / 3.0 );
    
     // Height Temp 4
    EstimateTempHeight( i_dt );
    // Vstar3, Astar3
    EstimateVelStar( i_dt );
    A_from_H( StateHeightTmp, StateVelStar, i_dt );
    // Accumulate
    AccumulateEstimate( i_dt / 6.0 );
    
    // Final boundary conditions on height and vel
    EnforceBoundaryConditions( StateHeight );
    EnforceBoundaryConditions( StateVel );

    // Update current time.
    StateCurrentTime += i_dt;
}

void DrawState()
{
    float OffsetY = 0.5 * ( float )WindowHeight;
    for ( int i = 0; i < ArraySize; ++i )
    {
        float SimX = DX *  ( float )i;
        float PixelsX = ( float )( i * PixelsPerCell );
        float SimY = State[StateHeight][i];
        float PixelsY = SimY * (( float )PixelsPerCell ) / DX;
        float PixelsMinY = OffsetY - PixelsY;
        float PixelsHeight = (( float )WindowHeight) - PixelsMinY;

        fill( 0.0, 0.0, 1.0 );   
        rect( PixelsX, OffsetY - PixelsY, PixelsPerCell, PixelsHeight );
    }
}

void draw()
{
    background( 0.9 );
   
    GetInput();
   
    TimeStep( 1.0 / 24.0 );

    DrawState();
    
    // Label.
    fill( 1.0 );
    text( "Wave Equation RK4 - With Input", 10, 30 );
}

// Reset function. If the key 'r' is released in the display, 
// copy the initial state to the state.
void keyReleased()
{
    if ( key == 114 )
    {
        SetInitialState();
    }  
}
