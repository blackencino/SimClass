//-*****************************************************************************
// Copyright (c) 2011-2013 Christopher Jon Horvath. All rights reserved.
//-*****************************************************************************

//-*****************************************************************************
// GLOBAL VARIABLES
//-*****************************************************************************

// Grid resolution per side. Rectangular
int NX = 128;
int NY = 128;

// The size of each grid cell, in pixels.
// This is for drawing
int CellSize = 5;

// The rate at which we inject density
// into the grid by painting with the
// mouse.
float EmissionRate = .01;

// The rate at which density
// diffuses (dissipates)
float D_viscosity = 0.000001;

// The rate at which velocity
// dissipates  
float V_viscosity = 0.00001;

// Our time step 
float DT = 1; 

// A scale on input velocity 
float VScale = 0.1 / ( float )( NX * CellSize );

// Our Window will be made of "gridRes" cells,
// where each cell is "cellSize" pixels big.
int WindowWidth = NX * CellSize;
int WindowHeight = NY * CellSize;

// Our simulation grids (Our State) will be one cell larger in each
// dimension to accomodate boundary conditions.
int GX = NX+2;
int GY = NY+2;

// The length of all of our (one-dimensional)
// arrays. We use 1d arrays rather than matrices
// mostly for efficiency reasons.
int GridArraySize = GX*GY;

// Here are our Velocity Arrays. 
// We use the standard convention of calling 
// X velocity "U" and  Y velocity "V".
// We have two sets of arrays - one which is our "current" array,
// and one which represents the previous array. We'll be clever with our
// indexing to support this.

float[][] GridsU = new float[2][GridArraySize];
float[][] GridsV = new float[2][GridArraySize];

// Density arrays
float[][] GridsDen = new float[2][GridArraySize];

// input velocities
float[] GridsInputU = new float[GridArraySize];
float[] GridsInputV = new float[GridArraySize];

// input source (or 'density')
float[] GridsInputDensity = new float[GridArraySize];

float v_strokeAlpha = 0.5;

// INDEX OF CURRENT STATE
int PREV = 0;
int CURR = 1;

//-*****************************************************************************
// Initialize everything
//-*****************************************************************************
void setup() 
{
    noSmooth();
    // Set up normalized colors.
    colorMode(RGB, 1.0);
    
    // The scale of our window
    size( WindowWidth, WindowHeight );
    
    // init density. Doing the full init.
    for ( int j = 0; j < GY; j++ )
    {
        float fj = j;
        for ( int i = 0; i < GX; i++ )
        {
            float fi = i;
            GridsDensity[CURR][IX(i,j)] = noise( fi / 16.0, fj / 16.0 );
        }
    }
}

//-*****************************************************************************
// Time Step
//-*****************************************************************************
void FluidTimeStep()
{
    // Put current into previous.
    SwapFluids();
    SwapDensities();

    // Get External Input

    // Solve Velocities
    IntegrateExternalForces();
    IntegrateViscosity();
    EnforceIncompressibility();
    AdvectFluid();
    EnforceIncompressibility();

    // Solve Densities
    DiffuseDensities();
    IntegrateExternalDensities();
    AdvectDensities();
}

//-*****************************************************************************
// Draw Everything
//-*****************************************************************************
void draw() 
{
    background( 1 );
    
    FluidTimeStep();

    drawDensity();
}

//-*****************************************************************************
// Draw Scalar Field
//-*****************************************************************************
void DrawScalarField( float[][] i_field, int I )
{
    loadPixels();
    for ( int wj = 0; wj < height; ++wj )
    {
        int j = ( int )( wj / CellSize );
        if ( j < 0 ) { j = 0; }
        else if ( j >= GY ) { j = GY-1; }

        for ( int wi = 0; wi < width; ++wi )
        {
            int i = ( int )( wi / CellSize );
            if ( i < 0 ) { i = 0; }
            else if ( i >= GX ) { i = GX-1; }

            float d = i_field[I][IX(i,j)];

            //colorMode(HSB, 1);
            //float h = map(d, 0, 1, 0.675, .55);
            //float s = map(d, 0, 1, 0, 1);
            //float b = map(d, 0, 1, 1, .25);

            pixels[ wi + wj*width ] = color( d, 1-d, 0.0 );
        }
    }
    updatePixels();
}

//-*****************************************************************************
// Draw Velocity Field
//-*****************************************************************************
void DrawVelocityField( float[][] i_velU, float[][] i_velV, int I )
{
    if ( v_strokeAlpha > 0.001 ) 
    {
        for ( int j = 1; j <= NY; ++j ) 
        {
            for ( int i = 1; i <= NX; ++i ) 
            {
                float cellMidpointX = i + 0.5;
                float cellMidpointY = j + 0.5;

                cellMidpointX *= cellSize;
                cellMidpointY *= cellSize;

                float v_x = cellMidpointX + input_u[IX(i,j)]/vScale;
                float v_y = cellMidpointY + input_v[IX(i,j)]/vScale;

                float dist = pow((v_x - cellMidpointX)/cellSize,2.0) + 
                 pow((v_y - cellMidpointY)/cellSize, 2.0);

                float vmag = map(dist, 0, 5, 0.1, 1);


                colorMode(HSB, 1);
                float h = map(vmag, 0, 1, 0, .025);
                float s = map(vmag, 0, 1, 1, .9);
                float b = map(vmag, 0, 1, 0, 1);

                stroke(h, s, b, v_strokeAlpha);

                line( cellMidpointX, cellMidpointY, v_x, v_y );
                colorMode(RGB, 1);


            } 
        }

        v_strokeAlpha *= 0.98;
    }
}


// --------------------------------
// Draw the velocity field
// --------------------------------
void drawVelocity()
{
  
  if (v_strokeAlpha > 0.00001) {
    for ( int i = 0 ;i < n; i++ ) {
      for ( int j = 0; j < n; j++ ) {
      
        float cellMidpointX = i + 0.5;
        float cellMidpointY = j + 0.5;
      
        cellMidpointX *= cellSize;
        cellMidpointY *= cellSize;
      
        float v_x = cellMidpointX + u[IX(i,j)]/vScale;
        float v_y = cellMidpointY + v[IX(i,j)]/vScale;
      
        float dist = pow((v_x - cellMidpointX)/cellSize,2.0) + 
                     pow((v_y - cellMidpointY)/cellSize, 2.0);
      
        float vmag = map(dist, 0, 50, 0.1, 1);
        
        colorMode(HSB, 1);
        float h = map(vmag, 0, 1, 0, .025);
        float s = map(vmag, 0, 1, 1, .9);
        float b = map(vmag, 0, 1, 0, .5);

        stroke(h, s, b, v_strokeAlpha);

        line( cellMidpointX, cellMidpointY, v_x, v_y );
        colorMode(RGB, 1);
        
      
      } 
    }
  }
}



// --------------------------------
// Add density based on mouse clicking
// --------------------------------
void GetInputSourceDensity()
{
    if ( mousePressed == true && mouseButton == LEFT ) 
    {
        for ( int j = 0; j < GY; ++j ) 
        {
            float cellMidPointY = CellSize * ( 0.5 + ( float )j );

            for ( int i = 0; i < GX; ++i ) 
            {
                float cellMidPointX = CellSize * ( 0.5 * ( float )i );

                float r = dist( mouseX, mouseY, cellMidPointX, cellMidPointY );
                float er2 = sq( r / EmissionRadius );

                GridsInputDensity[0][IX(i,j)] += 
                    emissionRate * constrain(1.0/dist, 0, 1);
            } 
        }
    } 
    else 
    {
        for ( int a = 0; a < ArraySize; ++a ) 
        {
            GridsInputDensity[0][a] = 0.0;
        }
    } 
}

// --------------------------------
// Derive the velocity of the mouse
// when it's clicked; add the resulting
// velocities to the grid velocity.
// --------------------------------
void getInputSourceVelocity ()
{

  if (mousePressed == true && mouseButton == RIGHT) {
    
    v_strokeAlpha = 0.5;  
    
    float m_vx = mouseX - pmouseX;
    float m_vy = mouseY - pmouseY;
  
  
    for ( int i = 0 ;i < n; i++ ) {
      for ( int j = 0; j < n; j++ ) {
        
        float cellMidpointX = i + 0.5;
        float cellMidpointY = j + 0.5;
        
        cellMidpointX *= cellSize;
        cellMidpointY *= cellSize;
      
        float dist = pow((mouseX - cellMidpointX)/cellSize, 2.0) + 
                     pow((mouseY - cellMidpointY)/cellSize, 2.0);
                     
        input_u[IX(i,j)] += m_vx * constrain(1.0/dist, 0, 1) * vScale;
        input_v[IX(i,j)] += m_vy * constrain(1.0/dist, 0, 1) * vScale;
        
      } 
    }
  } else {
    for ( int i = 0 ;i < size; i++ ) {
        input_u[i] *= 0.9;
        input_v[i] *= 0.9;
    }
  }
    
}


// -----------------------------------------------------------------
// Fluid Solver Stuff

// --------------------------------
// Add input density
// --------------------------------
void cool_density()
{
  for ( int i = 0 ;i < size; i++ ) {
      den[i] *= 0.999;
  }
}


// --------------------------------
// Add input density
// --------------------------------
void addSource(float[] x, float[] s)
{
  for ( int i = 0 ;i < size; i++ ) {
      x[i] += dt * s[i];
  }
}


// --------------------------------
// Diffuse (dissipate) density
// --------------------------------
void diffuse(float[] x, float[] x0, float visc, int b)
{
  int i, j, k;
  float a = dt * visc * n * n;
  
  
  for ( k=0 ; k<20 ; k++ ) {
    for ( i=1 ; i<=n ; i++ ) {
      for ( j=1 ; j<=n ; j++ ) {
        x[IX(i,j)] = (x0[IX(i,j)] + 
                        a * (x[IX(i-1,j)] + 
                             x[IX(i+1,j)] + 
                             x[IX(i,j-1)] + 
                             x[IX(i,j+1)] ) 
                        ) / ( 1 + 4*a);
      }
    }
    set_boundary( b, x );
  }
}

// --------------------------------
// Advect density by velocity
// --------------------------------
void advect(float[] d, float[] d0, float[] inU, float[] inV, int b)
{
  int i, j, i0, j0, i1, j1;
  float x, y, s0, t0, s1, t1, dt0;
  dt0 = dt * n;

  for ( i=1; i<=n; i++ ) {
    for ( j=1; j<=n; j++ ) {
    
      x = i - dt0 * inU[IX(i,j)];
      y = j - dt0 * inV[IX(i,j)];
      
      if (x < 0.5) x = 0.5; 
      if (x > n + 0.5) x = n + 0.5;
      i0=(int)x; 
      i1=i0+1;

      if ( y < 0.5) y = 0.5; 
      if ( y > n + 0.5) y = n + 0.5; 
      j0=(int)y; 
      j1=j0+1;
      
      s1 = x-i0;
      s0 = 1-s1; 
      t1 = y-j0; 
      t0 = 1-t1;
      
      d[IX(i,j)] = s0 * ( t0 * d0[IX(i0,j0)] + 
                          t1 * d0[IX(i0,j1)]) +
                   s1 * ( t0 * d0[IX(i1,j0)] + 
                          t1 * d0[IX(i1,j1)]) ;
    }
  }
  
  //set_boundary( b, d );

}

//-*****************************************************************************
// Compute Divergence.
// 
// Compute the divergence. This is the amount of mass entering
// or exiting each cell. In an incompressible fluid, divergence is
// zero - and that's what we're solving for!
//
// It is equal to the partial derivative of the x-velocity in the x direction
// (dU/dx),
// plus the partial derivative of the y-velocity in the y direction
// (dV/dy).
// We're using a central differencing scheme for computing the derivatives 
// here.
//-*****************************************************************************
void ComputeDivergence()
{
    for ( int j = 1; j <= NY; ++j )
    {
        for ( int i = 1; i <= NX; ++i )
        {
            float twoDU = U[CURR][IX(i+1,j)] - U[CURR][IX(i-1,j)];
            float twoDV = V[CURR][IX(i,j+1)] - V[CURR][IX(i,j-1)];
            GridsDivergence[CURR][IX(i,j)] = 
                ( twoDU / TwoDXY ) + ( twoDV / TwoDXY );
        }
    }
    
    // Compute Divergence Boundary conditions.
    EnforceBoundaryConditions( BC_NoNegate, GridsDivergence, CURR );
}

//-*****************************************************************************
// We're using the Projection Method for solving an incompressible flow.
//
// Information about this method can be found here:
// http://en.wikipedia.org/wiki/Projection_method_(fluid_dynamics)
//
// But, essentially, what it means is that we're breaking the solving
// of the equation up into separate steps, roughly:
// Add External Forces
// Advect
// Compute Divergence
// Find Pressure Which Will Correct Divergence
// Adjust Velocity By Negative Gradient Of Pressure.
//
// This function implements the Find Pressure part, by iteratively
// solving for a pressure which satisifies (at every point in the grid)
// the relationship:  Laplacian( Pressure ) = Divergence.
// The Laplacian is a differential operator which is equal to the
// second partial derivative with respect to x plus the second partial
// derivative with respect to y.  This can be discretized in the following
// way:
//
// SecondDerivX( P, i, j ) = ( ( P[i+1,j] - P[i,j] ) / DXY ) -
//                             ( P[i,j] - P[i-1,j] ) / DXY ) ) / DXY
// SecondDerivY( P, i, j ) = ( ( P[i,j+1] - P[i,j] ) / DXY ) -
//                             ( P[i,j] - P[i,j-1] ) / DXY ) ) / DXY
// 
// Laplacian( P, i, j ) = SecondDerivX( P, i, j ) + SecondDerivY( P, i, j )
// 
// This simplifies to:
// (( P[i,j-1] + P[i-1,j] + P[i+1,j] + P[i,j+1] ) - 4 P[i,j] ) / DXY^2
// 
// The idea with jacobi iteration is to assume, at every i,j point, that
// all the other values are constant, and just solve for one's self.
// This needs to be done into a new destination array, otherwise you'll be
// writing over your data as you're computing it. There are simulation
// techniques such as red-black gauss-seidel which can compute the data
// in-place, but we'll just swap between an old and a new pressure.
//
// Given the above expression for the lapacian at i,j, and the relationship
// Laplacian( P ) = Divergence
// We set the above expression to Divergence[i,j], and then rearrange the
// equation so that we're solving for P[i,j]:
//
// (( Pdown + Pleft + Pright + Pup ) - 4 Pcen )/H^2 = Div
// ( Pdown + Pleft + Pright + Pup ) - 4 Pcen  = H^2 * Div
// -4 Pcen = ( H^2 * Div ) - ( Pdown + Pleft + Pright + Pup )
// 4 Pcen = ( Pdown + Pleft + Pright + Pup ) - ( H^2 * Div )
// Pcen = ( ( Pdown + Pleft + Pright + Pup ) - ( H^2 * Div ) ) / 4
//
// P[i,j] = -( DXY^2 * Divergence[i,j] + 
//          (( P[i,j-1] + P[i-1,j] + P[i+1,j] + P[i,j+1] )) ) / 4
// 
//-*****************************************************************************
void ComputePressureViaJacobiIterations()
{
    // Init array indices.
    PPREV = 0;
    PCURR = 1;

    // Init the pressure curr to zero. It will be swapped into the prev
    // location in the loop below.
    for ( int a = 0; a < ArraySize; ++a )
    {
        GridsPressure[PCURR][a] = 0.0;
    }
    
    // Iterate 20 times, improving the pressure current from the pressure
    // prev.
    for ( int iter = 0; iter < 20; ++iter )
    {
        // Swap the indices of the current & previous pressure arrays.
        int tmp = PCURR;
        PCURR = PPREV;
        PPREV = tmp;

        // Do a single jacobi iteration to compute the current pressure
        // from the previous pressure.
        for ( int j = 1; j <= NY; ++j )
        {
            for ( int i = 1; i <= NX; ++i )
            {
                GridsPressure[PCURR][IX(i,j)] =
                ( ( GridsPressure[PPREV][IX(i,j-1)] +
                    GridsPressure[PPREV][IX(i-1,j)] +
                    GridsPressure[PPREV][IX(i+1,j)] +
                    GridsPressure[PPREV][IX(i,j+1)] ) -
                  ( DXY * DXY * GridsDivergence[CURR][IX(i,j)] ) ) / 4.0;
            }
        }

        // Okay we've solved for PCURR. Enforce boundary conditions on it,
        // without negating in any direction.
        EnforceBoundaryConditions( BC_NoNegate, GridsPressure, PCURR );
    }
}

//-*****************************************************************************
// Apply Negative Gradient of Pressure to Velocity
void ApplyNegativeGradientOfPressureToVelocity()
{ 
    for ( int j = 1; j <= NY; ++j )
    {
        for ( int i = 1; i <= NX; ++i )
        {
            float twoDPx = GridsPressure[PCURR][IX(i+1,j)] - 
                           GridsPressure[PCURR][IX(i-1,j)];
            float twoDPy = GridsPressure[PCURR][IX(i,j+1)] - 
                           GridsPressure[PCURR][IX(i,j-1)];
        
            GridsU[CURR][IX(i,j)] -= twoDPx / TwoDXY;
            GridsV[CURR][IX(i,j)] -= twoDPy / TwoDXY;
        }
    }

    // And apply boundary conditions. The U velocities are negated horizonally,
    // and the V velocities are negated vertically. This makes the fluid
    // reflect off the boundaries.
    EnforceBoundaryConditions( BC_NegateX, GridsU, CURR );
    EnforceBoundaryConditions( BC_NegateY, GridsV, CURR );
}

//-*****************************************************************************
// The boundary conditions are enforced on the I selector of the arrays.
// There are three types of boundary condition application - no negation,
// just copying at the boundary, then negating in the x-direction only,
// then negating in the y-direction only.
void EnforceBoundaryConditions( int i_bType, float[][] io_array, int I )
{
    // Copy the bottom row from the one above it. If the boundary type
    // is '2', that means negate the values.
    // Do the same for the topmost row and the one beneath it.
    for ( int i = 1; i <= NX; ++i )
    {
        io_array[I][IX( i, 0  )] = 
            i_btype==BC_NegateY ? -io_array[I][IX( i, 1 )];
        io_array[I][IX( i, NY+1 )] = 
            i_btype==BC_NegateY ? -io_array[I][IX( i, NY )];
    }

    // Copy the left col from the one to the right of it. If the boundary type
    // is '1', that means negate the values.
    // Do the same for the rightmost col and the one to the left of it.
    for ( int j = 1; j <= NY; ++i )
    {
        io_array[I][IX( 0, j  )]   = 
            i_btype==BC_NegateX ? -io_array[I][IX( 1, j )];
        io_array[I][IX( NX+1, j )] = 
            i_btype==BC_NegateY ? -io_array[I][IX( NX, j )];
    }

    // Get each corner by averaging the two boundary values adjacent.
    io_array[I][IX(0,0)]     
        = 0.5 * ( io_array[I][IX(1,0)] + io_array[I][IX(0,1)] );
    io_array[I][IX(0,NY+1)]   
        = 0.5 * ( io_array[I][IX(1,NY+1)] + io_array[I][IX(0,NY)] );
    io_array[I][IX(NX+1,0)]   
        = 0.5 * ( io_array[I][IX(NX,0)] + io_array[I][IX(NX+1,1)] );
    io_array[I][IX(NX+1,NY+1)] 
        = 0.5 * ( io_array[I][IX(NX,NY+1)] + io_array[I][IX(NX+1,NY)]);
}

//-*****************************************************************************
// Utility functions
//-*****************************************************************************


// Index an element of our grid
int IX( int i, int j )
{
    return ( i + GX*j ); 
}

void SWAP(float[] x, float[] x0)
{

  float tmp = 0;
  for ( int i = 0 ;i < size; i++ ) {
      tmp = x0[i];
      x0[i] = x[i];
      x[i] = tmp;
  }
}

