float WaveSpeed = 0.5;

float WorldSize = 10.0;
int ArraySize = 128;
float DX = WorldSize / ArraySize;

float LX = WorldSize;
float LY = WorldSize / 2.0;

int StateSize = 4;
float[][] State = new float[StateSize][ArraySize];
int StateHeight = 0;
int StateHeightPrev = 1;
int StateHeightPrevPrev = 2;
int StateJacobiTmp = 3;

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

void EnforceBoundaryConditions( int io_a )
{
    State[io_a][0] = ( 2.0 * State[io_a][1] ) - State[io_a][2];
    State[io_a][ArraySize-1] = ( 2.0 * State[io_a][ArraySize-2] ) - State[io_a][ArraySize-3];
    //State[io_a][0] = State[io_a][1];
    //State[io_a][ArraySize-1] = State[io_a][ArraySize-2];
    
    //if ( InputActive && io_a == StateHeight )
    //{
    //    State[io_a][InputIndex] = InputHeight;
    //}
}

void CopyArray( int i_src, int o_dst )
{
    for ( int i = 0; i < ArraySize; ++i )
    {
        State[o_dst][i] = State[i_src][i];
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
    }
    EnforceBoundaryConditions( StateHeight );
    CopyArray( StateHeight, StateHeightPrev );
    CopyArray( StateHeight, StateHeightPrevPrev );
    InputActive = false;
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

void RotateState()
{
    int tmp = StateHeightPrevPrev;
    StateHeightPrevPrev = StateHeightPrev;
    StateHeightPrev = StateHeight;
    StateHeight = tmp;
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

// Jacobi iteration to get temp height.
void JacobiIterationHeight( int i_hOld, int o_hNew, float i_dt )
{
    float aij = -sq( WaveSpeed ) * sq( i_dt ) / sq( DX );
    float aii = 1.0 + ( 2.0 * sq( WaveSpeed ) * sq( i_dt ) / sq( DX ) );
    
    for ( int i = 1; i < ArraySize-1; ++i )
    {
        float hLeft = State[i_hOld][i-1];
        float hRight = State[i_hOld][i+1];
        
        float hp = State[StateHeightPrev][i];
        float hpp = State[StateHeightPrevPrev][i];
        
        float bi = ( 2.0 * hp ) - hpp;
        
        float sumOffDiag = ( aij * hLeft ) + ( aij * hRight );
        
        State[o_hNew][i] = ( bi - sumOffDiag ) / aii;
    }
    
    EnforceBoundaryConditions( o_hNew );
}

void JacobiSolveHeight( float i_dt )
{  
    // First iteration solve from StateHeightPrev into StateJacobiTmp
    JacobiIterationHeight( StateHeightPrev, StateHeight, i_dt );
    
    // All other iterations solve from StateJacobiTmp into StateHeight
    for ( int iter = 0; iter < 20; ++iter )
    {
        int tmp = StateHeight;
        StateHeight = StateJacobiTmp;
        StateJacobiTmp = tmp;
        
        JacobiIterationHeight( StateJacobiTmp, StateHeight, i_dt );
    }
}

// Time Step function.
void TimeStep( float i_dt )
{
    // Rotate state.
    RotateState();
    
    // Jacobi-solve height.
    JacobiSolveHeight( i_dt );
    
    if ( InputActive )
    {
        State[StateHeight][InputIndex] = InputHeight;
    }
    
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
