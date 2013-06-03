float WaveSpeed = 0.5;

float WorldSize = 10.0;
int ArraySize = 128;
float DX = WorldSize / ArraySize;

float LX = WorldSize;
float LY = WorldSize / 2.0;

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

int PixelsPerCell = 4;

int WindowWidth = PixelsPerCell * ArraySize;
int WindowHeight = WindowWidth / 2;

boolean InputActive = false;
int InputIndex = 0;
float InputHeight = 0;

float snoise( float x )
{
   return ( 2.0 * noise( x )) - 1.0;
}

float anoise( float x )
{
   return ( -2.0 * abs( snoise( x )) ) + 1.0;
}

void EnforceAccelBoundaryConditions( int io_a )
{
    State[io_a][0] = 0.0;
    State[io_a][ArraySize-1] = 0.0;
}

void EnforceVelBoundaryConditions( int io_v )
{
    State[io_v][0] = State[io_v][1];
    State[io_v][ArraySize-1] = State[io_v][ArraySize-2];
}

void EnforceHeightBoundaryConditions( int io_h )
{
    //State[io_h][0] = ( 2.0 * State[io_h][1] ) - State[io_h][2];
    //State[io_h][ArraySize-1] = ( 2.0 * State[io_h][ArraySize-2] ) - State[io_h][ArraySize-3];
    State[io_h][0] = State[io_h][1];
    State[io_h][ArraySize-1] = State[io_h][ArraySize-2];
    
    if ( InputActive )
    {
        State[io_h][InputIndex] = InputHeight;
    }
}

void CopyArray( int i_src, int o_dst )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[o_dst][i] = State[i_src][i];
    }
}

void FillArray( int o_a, float i_val )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[o_a][i] = i_val;
    }
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
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceVelBoundaryConditions( StateVel );
    StateCurrentTime = 0.0;
    
    CopyArray( StateHeight, StateHeightPrev );
    CopyArray( StateVel, StateVelPrev );
    InputActive = false;
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

void GetInput()
{
    InputActive = false;
    
    if ( mousePressed && mouseButton == LEFT )
    {
        float mouseCellX = mouseX / ( ( float )PixelsPerCell );
        float mouseCellY = ( (height/2) - mouseY ) / ( ( float )PixelsPerCell );
        float simY = mouseCellY * DX;
        
        int iX = ( int )floor( mouseCellX + 0.5 );
        if ( iX > 0 && iX < ArraySize-1 )
        {
            InputIndex = iX;
            InputHeight = simY;
            InputActive = true;
        }
    }
}

// Estimate height star
void EstimateHeightStar( float i_dt )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[StateHeightStar][i] = State[StateHeightPrev][i] + 
                ( i_dt * State[StateVelStar][i] );
    }
    EnforceHeightBoundaryConditions( StateHeightStar );
}

// Estimate vel star
void EstimateVelStar( float i_dt )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[StateVelStar][i] = State[StateVelPrev][i] + 
                ( i_dt * State[StateAccelStar][i] );
    }
    EnforceVelBoundaryConditions( StateVelStar );
}

// Jacobi iteration to get temp acceleration
void JacobiIterationAccel( int i_aOld, int o_aNew, int i_hStar, float i_dt )
{
    float aij = -sq( WaveSpeed ) * sq( i_dt ) / sq( DX );
    float aii = 1.0 + ( 2.0 * sq( WaveSpeed ) * sq( i_dt ) / sq( DX ) );
    float bgain = sq( WaveSpeed ) / sq( DX );
    
    for ( int i = 1; i < ArraySize-1; ++i )
    {
        float aLeft = State[i_aOld][i-1];
        float aRight = State[i_aOld][i+1];
        
        float hStarLeft = State[i_hStar][i-1];
        float hStarCen = State[i_hStar][i];
        float hStarRight = State[i_hStar][i+1];
        
        float bi = bgain * ( hStarLeft + hStarRight - ( 2.0 * hStarCen ) );
        
        float sumOffDiag = ( aij * aLeft ) + ( aij * aRight );
        
        State[o_aNew][i] = ( bi - sumOffDiag ) / aii;
    }
    
    EnforceAccelBoundaryConditions( o_aNew );
}

// Solve for acceleration.
void JacobiSolveAccel( int i_hStar, float i_dt )
{  
    // Initialize acceleration to zero.
    FillArray( StateAccelStar, 0.0 );
    
    // Solve from StateJacobiTmp into StateAccel
    for ( int iter = 0; iter < 20; ++iter )
    {
        int tmp = StateAccelStar;
        StateAccelStar = StateJacobiTmp;
        StateJacobiTmp = tmp;
        
        JacobiIterationAccel( StateJacobiTmp, StateAccelStar, i_hStar, i_dt );
    }
}

void EstimateAccelStar( float i_dt )
{
    JacobiSolveAccel( StateHeightStar, i_dt );
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

// Time Step function.
void TimeStep( float i_dt )
{
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
    EnforceVelBoundaryConditions( StateVel );

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
    text( "Wave Equation Stable Jacobi RK4 - With Input", 10, 30 );
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
