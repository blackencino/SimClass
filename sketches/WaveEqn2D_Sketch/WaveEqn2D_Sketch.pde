float WaveSpeed = 0.5;

float WorldSize = 10.0;
int NX = 64;
int NY = 64;
int ArraySize = NX * NY;
float DX = WorldSize / NX;
float DY = WorldSize / NY;

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

PImage StateImage = createImage( NX, NY, RGB );


// Index an element of a grid in the state array
int IX( int i, int j )
{
    return ( i + NX*j ); 
}

float snoise( float x, float y )
{
   return ( 2.0 * noise( x, y )) - 1.0;
}

float anoise( float x, float y )
{
   return ( -2.0 * abs( snoise( x, y )) ) + 1.0;
}

void EnforceAccelBoundaryConditions( int io_a )
{
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

void EnforceVelBoundaryConditions( int io_v )
{
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

void EnforceHeightBoundaryConditions( int io_h )
{
    for (int j = 0; j < NY; ++j) {
        if (j == 0) {
            for (int i = 0; i < NX; ++i) {
                State[io_h][IX(i,0)] = State[io_h][IX(i,1)];
            }
        } else if (j == (NY-1)) {    
            for (int i = 0; i < NX; ++i) {
                State[io_h][IX(i,NY-1)] = State[io_h][IX(i,NY-2)];
            }
        }

        State[io_h][IX(0,j)] = State[io_h][IX(1,j)];
        State[io_h][IX(NX-1,j)] = State[io_h][IX(NX-2,j)];
    }
    
    if ( InputActive )
    {
        float center_weight = 1.0;
        float edge_weight = 0.25;
        float corner_weight = edge_weight * edge_weight;
        
        State[io_h][IX(InputIndexX-1, InputIndexY-1)] = 
          lerp(State[io_h][IX(InputIndexX-1, InputIndexY-1)], 1.0, corner_weight);        
        State[io_h][IX(InputIndexX, InputIndexY-1)] = 
          lerp(State[io_h][IX(InputIndexX, InputIndexY-1)], 1.0, edge_weight);
        State[io_h][IX(InputIndexX+1, InputIndexY-1)] = 
          lerp(State[io_h][IX(InputIndexX+1, InputIndexY-1)], 1.0, corner_weight);
        
        State[io_h][IX(InputIndexX-1, InputIndexY)] = 
          lerp(State[io_h][IX(InputIndexX-1, InputIndexY)], 1.0, edge_weight);        
        State[io_h][IX(InputIndexX, InputIndexY)] = 
          lerp(State[io_h][IX(InputIndexX, InputIndexY)], 1.0, center_weight);
        State[io_h][IX(InputIndexX+1, InputIndexY)] = 
          lerp(State[io_h][IX(InputIndexX+1, InputIndexY)], 1.0, edge_weight);
        
        State[io_h][IX(InputIndexX-1, InputIndexY+1)] = 
          lerp(State[io_h][IX(InputIndexX-1, InputIndexY+1)], 1.0, corner_weight);        
        State[io_h][IX(InputIndexX, InputIndexY+1)] = 
          lerp(State[io_h][IX(InputIndexX, InputIndexY+1)], 1.0, edge_weight);
        State[io_h][IX(InputIndexX+1, InputIndexY+1)] = 
          lerp(State[io_h][IX(InputIndexX+1, InputIndexY+1)], 1.0, corner_weight);
        
        //State[io_h][IX(InputIndexX, InputIndexY)] = InputHeight;
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
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            float worldX = 2341.17 + DX * ( float )i;
            float worldY = 9911.44 + DY * ( float )j;

            State[StateHeight][IX(i,j)] = 
                0.5 * anoise( worldX * 0.0625, worldY * 0.0625 ) +
                0.4 * anoise( worldX * 0.125, worldY * 0.125 ) +
                0.3 * anoise( worldX * 0.25, worldY * 0.25 ) +
                0.2 * anoise( worldX * 0.5, worldY * 0.5 );
            State[StateVel][IX(i,j)] = 0.0;
        }
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
        float mouseCellY = mouseY / ( ( float )PixelsPerCell );
        
        int iX = ( int )floor( mouseCellX + 0.5 );
        int iY = ( int )floor( mouseCellY + 0.5 );
        if ( iX > 0 && iX < NX-1 &&
             iY > 0 && iY < NY-1 )
        {
            InputIndexX = iX;
            InputIndexY = iY;
            InputHeight = 1.0;
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
    float kappa = sq( WaveSpeed ) * sq( i_dt ) / sq( DX );
    float gamma = sq( WaveSpeed ) / sq( DX );
    float coefficient_diagonal = 1.0 + 4.0 * kappa;
    float coefficient_off_diagonal = -kappa;
    float rhs_coefficient_gain = gamma;

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

            float rhs_coefficient = rhs_coefficient_gain * 
                (h_star_left + h_star_right + h_star_down + h_star_up -
                 (4.0 * h_star_cen));

            float sum_off_diagonal = coefficient_off_diagonal * 
                (a_left + a_right + a_down + a_up);

            State[o_aNew][IX(i,j)] = (rhs_coefficient - sum_off_diagonal) /
                    coefficient_diagonal;
        }
    }
    
    EnforceAccelBoundaryConditions( o_aNew );
}

// Solve for acceleration.
void JacobiSolveAccel( int i_hStar, float i_dt )
{  
    // Initialize acceleration to zero.
    FillArray( StateAccelStar, 0.0 );
    
    // Solve from StateJacobiTmp into StateAccel
    for ( int iter = 0; iter < 10; ++iter )
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
void TimeStepFirstOrder( float i_dt )
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
    AccumulateEstimate( i_dt );

    // Final boundary conditions on height and vel
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceVelBoundaryConditions( StateVel );

    // Update current time.
    StateCurrentTime += i_dt;
}

// Time Step function.
void TimeStepRK2( float i_dt )
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
    AccumulateEstimate( i_dt / 2.0 );
    
    // 2
    EstimateVelStar( i_dt );
    EstimateHeightStar( i_dt );
    EstimateAccelStar( i_dt );
    AccumulateEstimate( i_dt / 2.0 );
    
    // Final boundary conditions on height and vel
    EnforceHeightBoundaryConditions( StateHeight );
    EnforceVelBoundaryConditions( StateVel );

    // Update current time.
    StateCurrentTime += i_dt;
}

// Time Step function.
void TimeStepRK4( float i_dt )
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

float schlick(float n2, float cos_theta) {
    float n1 = 1.0;
    float R0 = sq((n1 - n2)/(n1 + n2));
    return R0 + (1.0 - R0)*pow(1.0 - cos_theta, 5.0);  
}

float kSpecular(PVector In, PVector Nn, PVector Ln, float m )
{
    PVector Vn = PVector.mult(In, -1.0);
    PVector H = PVector.add(Ln, Vn);
    H.normalize();
    float d = Nn.dot(H); 
    d *= d;
    return pow( d, m/2 );
}

float kDiffuse(PVector Nn, PVector Ln)
{
    float d = Nn.dot(Ln);
    return constrain( d, 0, 1 );
}

//-*****************************************************************************
// DRAW HEIGHT FIELD
//-*****************************************************************************
void DrawHeightField( int i_field ) {
    float pixr, pixg, pixb;
    float d;
    StateImage.loadPixels();
    PVector dh_di = new PVector(0, 0, 0);
    PVector dh_dj = new PVector(0, 0, 0);
    PVector n = new PVector(0, 0, 0);
    PVector toSun = new PVector(1, 1, 2.5);
    toSun.normalize();
    PVector toMoon = new PVector(-1, -.5, 1.25);
    toMoon.normalize();
    PVector In = new PVector(0, 0, -1.0);
    float specular_intensity_r = 30.0;
    float specular_intensity_g = 30.0;
    float specular_intensity_b = 30.0;
    float diffuse_intensity_r = 0.05;
    float diffuse_intensity_g = 0.125;
    float diffuse_intensity_b = 0.0625;
    float sun_r = 1.0;
    float sun_g = 1.0;
    float sun_b = 1.0;
    float moon_r = 0.1;
    float moon_g = 0.1;
    float moon_b = 0.25;
    for (int j = 0; j < NY; ++j) {
        int j_neg = constrain(j-1, 0, NY-1);
        int j_pos = constrain(j+1, 0, NY-1);
        float dy;
        if (j == 0 || j == (NY-1)) {
            dy = DY;
        } else {
            dy = 2.0 * DY;
        }
        for (int i = 0; i < NX; ++i) {
            int i_neg = constrain(i-1, 0, NX-1);
            int i_pos = constrain(i+1, 0, NX-1);
            float dx;
            if (i == 0 || i == (NX-1)) {
                dx = DX;
            } else {
                dx = 2.0 * DX;
            }
            float h_left = State[i_field][IX(i_neg,j)];
            float h_right = State[i_field][IX(i_pos,j)];
            float h_down = State[i_field][IX(i,j_neg)];
            float h_up = State[i_field][IX(i,j_pos)];
            float h_cen = State[i_field][IX(i,j)];
            
            dh_di.set(dx, 0.0, h_right - h_left);
            dh_dj.set(0.0f, dy, h_up - h_down);
            n.set(dh_di.cross(dh_dj));
            n.normalize();
            
            float fresnel_d = schlick(1.31, n.z);
            float sun_spec = fresnel_d * kSpecular(In, n, toSun, 80);//sq(n.dot(toSun));
            float sun_spec_r = specular_intensity_r * sun_spec * sun_r;
            float sun_spec_g = specular_intensity_g * sun_spec * sun_g;
            float sun_spec_b = specular_intensity_b * sun_spec * sun_b;
            
            float sun_diff = kDiffuse(n, toSun);
            float sun_diff_r = diffuse_intensity_r * sun_diff * sun_r;
            float sun_diff_g = diffuse_intensity_g * sun_diff * sun_g;
            float sun_diff_b = diffuse_intensity_b * sun_diff * sun_b;
            
            
            float moon_spec = fresnel_d * kSpecular(In, n, toMoon, 50);//sq(n.dot(toSun));
            float moon_spec_r = specular_intensity_r * moon_spec * moon_r;
            float moon_spec_g = specular_intensity_g * moon_spec * moon_g;
            float moon_spec_b = specular_intensity_b * moon_spec * moon_b;
            
            float moon_diff = kDiffuse(n, toMoon);
            float moon_diff_r = diffuse_intensity_r * moon_diff * moon_r;
            float moon_diff_g = diffuse_intensity_g * moon_diff * moon_g;
            float moon_diff_b = diffuse_intensity_b * moon_diff * moon_b;
          
            //d = constrain(State[i_field][IX(i,j)], -1.0, 1.0);
            //d = 0.5 + (0.5 * d);

            //colorMode(HSB, 1);
            //float h = map(d, 0, 1, 0.675, .55);
            //float s = map(d, 0, 1, 0, 1);
            //float b = map(d, 0, 1, 1, .25);
            //pixr = 0.9 * ( 1.0 - d );
            //pixg = 0.9 * ( 1.0 - (d*d) );
            //pixb = 0.9 * ( 1.0 - (d*d*d) );
            pixr = (sun_spec_r + sun_diff_r + moon_spec_r + moon_diff_r);
            pixg = (sun_spec_g + sun_diff_g + moon_spec_g + moon_diff_g);
            pixb = (sun_spec_b + sun_diff_b + moon_spec_b + moon_diff_b);

            StateImage.pixels[IX(i,j)] = color( pixr, pixg, pixb );
        }
    }
    StateImage.updatePixels();
    image( StateImage, 0, 0, width, height );
}


//-*****************************************************************************
void draw()
{
    background( 0.5 );

    GetInput();
    TimeStepRK2( 1.0 / 24.0 );
    DrawHeightField(StateHeight);

    // Label.
    fill( 1.0 );
    text( "2D Wave Equation Stable Jacobi RK4", 10, 30 );
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

