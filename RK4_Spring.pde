<script type="application/processing" data-processing-target="pjsSpringRK4">
float Stiffness = 5.0;
float BobMass = 0.5;
int StateSize = 3;
float[] InitState = new float[StateSize];
float[] State = new float[StateSize];
int StateCurrentTime = 0;
int StatePositionX = 1;
int StateVelocityX = 2;

int WindowWidthHeight = 300;
float WorldSize = 2.0;
float PixelsPerMeter;
float OriginPixelsX;
float OriginPixelsY;

void setup()
{
    // Create initial state.
    InitState[StateCurrentTime] = 0.0;
    InitState[StatePositionX] = 0.65;
    InitState[StateVelocityX] = 0.0;

    // Copy initial state to current state.
    // notice that this does not need to know what the meaning of the
    // state elements is, and would work regardless of the state's size.
    for ( int i = 0; i < StateSize; ++i )
    {
        State[i] = InitState[i];
    }

    // Set up normalized colors.
    colorMode( RGB, 1.0 );
    
    // Set up the stroke color and width.
    stroke( 0.0 );
    //strokeWeight( 0.01 );
    
    // Create the window size, set up the transformation variables.
    size( WindowWidthHeight, WindowWidthHeight );
    PixelsPerMeter = (( float )WindowWidthHeight ) / WorldSize;
    OriginPixelsX = 0.5 * ( float )WindowWidthHeight;
    OriginPixelsY = 0.5 * ( float )WindowWidthHeight;

    textSize( 24 );
}

// Draw our State, with the unfortunate units conversion.
void DrawState()
{
    // Compute end of arm.
    float SpringEndX = PixelsPerMeter * State[StatePositionX];

    // Compute the CORRECT position.
    float sqrtKoverM = sqrt( Stiffness / BobMass );
    float x0 = InitState[StatePositionX];
    float v0 = InitState[StateVelocityX];
    float t = State[StateCurrentTime];
    float CorrectPositionX = ( x0 * cos( sqrtKoverM * t ) ) +
        ( ( v0 / sqrtKoverM ) * sin( sqrtKoverM + t ) );
    
    // Compute draw pos for "correct"
    float CorrectEndX = PixelsPerMeter * CorrectPositionX;

    // Draw the spring.
    strokeWeight( 1.0 );
    line( 0.0, 0.0, SpringEndX, 0.0 );
          
    // Draw the spring pivot
    fill( 0.0 );
    ellipse( 0.0, 0.0, 
             PixelsPerMeter * 0.03, 
             PixelsPerMeter * 0.03 );
    
    // Draw the spring bob
    fill( 1.0, 0.0, 0.0 );
    ellipse( SpringEndX, 0.0, 
             PixelsPerMeter * 0.1, 
             PixelsPerMeter * 0.1 );

    // Draw the correct bob in blue
    fill( 0.0, 0.0, 1.0 );
    ellipse( CorrectEndX, -PixelsPerMeter * 0.25,
             PixelsPerMeter * 0.1,
             PixelsPerMeter * 0.1 );
}

// Acceleration from Position.
float A_from_X( float i_x )
{
    return -( Stiffness / BobMass ) * i_x;
}

// Time Step function.
void TimeStep( float i_dt )
{
    float vStar1 = State[StateVelocityX];
    float aStar1 = A_from_X( State[StatePositionX] );

    float vStar2 = State[StateVelocityX] + ( ( i_dt / 2.0 ) * aStar1 );
    float xTmp2 = State[StatePositionX] + ( ( i_dt / 2.0 ) * vStar1 );
    float aStar2 = A_from_X( xTmp2 );

    float vStar3 = State[StateVelocityX] + ( ( i_dt / 2.0 ) * aStar2 );
    float xTmp3 = State[StatePositionX] + ( ( i_dt / 2.0 ) * vStar2 );
    float aStar3 = A_from_X( xTmp3 );

    float vStar4 = State[StateVelocityX] + ( i_dt * aStar3 );
    float xTmp4 = State[StatePositionX] + ( i_dt * vStar3 );
    float aStar4 = A_from_X( xTmp4 );

    State[StatePositionX] += ( i_dt / 6.0 ) * 
        ( vStar1 + (2.0*vStar2) + (2.0*vStar3) + vStar4 );
    State[StateVelocityX] += ( i_dt / 6.0 ) * 
        ( aStar1 + (2.0*aStar2) + (2.0*aStar3) + aStar4 );

    // Update current time.
    State[StateCurrentTime] += i_dt;
}

// Processing Draw function, called every time the screen refreshes.
void draw()
{
    // Time Step.
    TimeStep( 1.0/24.0 );

    // Clear the display to a constant color.
    background( 0.75 );

    // Label.
    fill( 1.0 );
    text( "RK4", 10, 30 );

    pushMatrix();

    // Translate to the origin.
    translate( OriginPixelsX, OriginPixelsY );

    // Draw the simulation
    DrawState();
}

// Reset function. If the key 'r' is released in the display, 
// copy the initial state to the state.
void keyReleased()
{
    if ( key == 114 )
    {
        for ( int i = 0; i < StateSize; ++i )
        {
            State[i] = InitState[i];
        }
    }  
}
