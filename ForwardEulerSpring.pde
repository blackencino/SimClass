float Stiffness = 1.0;
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
}

// Draw our State, with the unfortunate units conversion.
void DrawState()
{
    // Compute end of arm.
    float SpringEndX = PixelsPerMeter * State[StatePositionX];

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
}

// Time Step function.
void TimeStep( float i_dt )
{
    // Compute acceleration from current position.
    float A = ( -Stiffness / BobMass ) * State[StatePositionX];

    // Update position based on current velocity.
    State[StatePositionX] += i_dt * State[StateVelocityX];

    // Update velocity based on acceleration.
    State[StateVelocityX] += i_dt * A;

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

    // Translate to the origin.
    translate( OriginPixelsX, OriginPixelsY );

    // Draw the simulation
    DrawState();
}
