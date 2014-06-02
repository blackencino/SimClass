float PendulumLength = 1.0;
float PendulumInitTheta = radians( 30.0 );
float StatePendulumTheta = radians( 30.0 );
float StatePendulumDthetaDt = 0.0;

float WorldSize = 2.0;

float PixelsPerMeter;

float OriginPixelsX;
float OriginPixelsY;

int WindowWidthHeight = 300;

void setup()
{
    // Set up normalized colors.
    colorMode( RGB, 1.0 );
    
    // Set up the stroke color and width.
    stroke( 0.0 );
    //strokeWeight( 0.01 );
    
    // Create the window size, set up the transformation variables.
    size( WindowWidthHeight, WindowWidthHeight );
    PixelsPerMeter = (( float )WindowWidthHeight ) / WorldSize;
    OriginPixelsX = 0.5 * ( float )WindowWidthHeight;
    OriginPixelsY = 0.25 * ( float )WindowWidthHeight;
}

// The DrawSim function assumes that the coordinate space is that of the
// simulation - namely, meters, with the pendulum pivot placed at the origin.
// Draw the arm, pivot, and bob!
// There is currently a bug in processing.js which requires us to do the
// pixels-per-meter scaling ourselves.
void DrawSim()
{
    // Compute end of arm.
    float ArmEndX = PixelsPerMeter * PendulumLength * sin( StatePendulumTheta );
    float ArmEndY = PixelsPerMeter * PendulumLength * cos( StatePendulumTheta );
  
    // Draw the pendulum arm.
    //line( 0.0, 0.0, ArmEndX, ArmEndY );
    strokeWeight( 1.0 );
    line( 0.0, 0.0, ArmEndX, ArmEndY );
          
    // Draw the pendulum pivot
    fill( 0.0 );
    ellipse( 0.0, 0.0, 
             PixelsPerMeter * 0.03, 
             PixelsPerMeter * 0.03 );
    
    // Draw the pendulum bob
    fill( 1.0, 0.0, 0.0 );
    ellipse( ArmEndX, ArmEndY, 
             PixelsPerMeter * 0.1, 
             PixelsPerMeter * 0.1 );
}

// The draw function creates a transformation matrix between pixel space
// and simulation space, in meters, and then calls the DrawSim function.
// Unfortunately, there is currently a bug in processing.js with concatenated
// transformation matrices, so we have to do the coordinate scaling ourselves
// in the draw function.
void draw()
{
    background( 0.75 );

    // Translate to the origin.
    translate( OriginPixelsX, OriginPixelsY );

    // Draw the simulation
    DrawSim();
}
