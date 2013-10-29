float WorldSize = 10.0;
int ArraySize = 128;
float DX = WorldSize / ArraySize;

float[] StateHeight = new float[ArraySize];

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

void setup()
{
    noiseSeed( 0 );
    for ( int i = 0; i < ArraySize; ++i )
    {
        float worldX = 2341.17 + DX * ( float )i;
        StateHeight[i] = 0.5 * anoise( worldX * 0.0625 ) +
                         0.4 * anoise( worldX * 0.125 ) +
                         0.3 * anoise( worldX * 0.25 ) +
                         0.2 * anoise( worldX * 0.5 );
    }
    
    size( WindowWidth, WindowHeight );
    
    colorMode( RGB, 1.0 );
    strokeWeight( 0.5 );
}

void DrawState()
{
    float OffsetY = 0.5 * ( float )WindowHeight;
    for ( int i = 0; i < ArraySize; ++i )
    {
        float SimX = DX *  ( float )i;
        float PixelsX = ( float )( i * PixelsPerCell );
        float SimY = StateHeight[i];
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
  
   DrawState();
}
