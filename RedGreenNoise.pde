// --------------------------------
// Init Global vars
// --------------------------------

  // --------------------------------
  // Some config options 
  // for our simulator
  //
  //Grid resolution per side
  int NX = 128;
  int NY = 128;
  int CellSize = 4;

  // Our Window will be made of "gridRes" cells,
  // where each cell is "cellSize" pixels big.
  int WindowWidth = CellSize * NX;
  int WindowHeight = CellSize * NY;

  // The length of all of our (one-dimensional)
  // arrays. We use 1d arrays rather than matrices
  // mostly for efficiency reasons.
  int ArraySize = NX * NY;
  
  // density arrays
  float[] Density = new float[ArraySize];

// --------------------------------
// Initialize everything
// --------------------------------
void setup() 
{
  //noSmooth();
  // Set up normalized colors.
  colorMode(RGB, 1.0);

  noiseSeed( 0 );
  
  // The scale of our window
  size( WindowWidth, WindowHeight );
  
  // init density
  for ( int j = 0; j < NY; j++ )
  {
      float fj = j;
      for ( int i = 0; i < NX; i++ )
      {
          float fi = i;
          Density[ i + (NX*j) ] = noise( fi / 16.0, fj / 16.0 );
      }
  }
}


// --------------------------------
// Draw the density field
// --------------------------------
void drawFloatArray( float[] i_values )
{
  loadPixels();
  for ( int wj = 0; wj < height; ++wj )
  {
      int j = ( int )( wj / CellSize );

      for ( int wi = 0; wi < width; ++wi )
      {
          int i = ( int )( wi / CellSize );

          float d = i_values[ i + (j*NX) ];
      
          //colorMode(HSB, 1);
          //float h = map(d, 0, 1, 0.675, .55);
          //float s = map(d, 0, 1, 0, 1);
          //float b = map(d, 0, 1, 1, .25);

          pixels[ wi + wj*width ] = color( d, 1-d, 0.0 );
      
    }
  }
  updatePixels();
}

// --------------------------------
// Draw Everything
// --------------------------------
void draw() {
  
  background(1);
  
  drawFloatArray( Density );

  stroke( 0 );
  fill( 0 );
  ellipse( mouseX, mouseY, 25, 25 );
}


