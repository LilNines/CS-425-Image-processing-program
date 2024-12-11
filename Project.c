// Term Project
// Xialei Fang
// 200434895
// -----------------------------
// 1. Apply Laplacian Operator to the “building.raw” and find all zero-crossing points.
// 2. Apply Sobel Operator to the “building.raw” image to obtain the gradient at all pixels.
// 3. Generate an edge image where a pixel is an edge point if it is a zero-crossing point and 
//    the gradient at this point is greater than or equal to a pre-specified threshold.
// 4. Implement the Hough transform algorithm to extract three longest linear structures 
//    from your edge image.
//
// command: 
//     gcc Project.c -o project -lm
//     ./project input_image.raw output_zero.raw output_gradient.raw output_edge_map.raw output_hough.raw output_final.raw
//      For Example: ./project building.raw zero.raw gradient.raw edge.raw hough.raw final.raw
// building.raw - 560 x 420

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ROWS 420   // This is hard-coded image resolution
#define COLS 560   // Dimension of the image
#define M 3        // For 3 × 3 Sobel operators

#define PI 3.14159265
#define rMax (int)(sqrt(ROWS * ROWS + COLS * COLS))
#define tMax 180

unsigned char in_buf[ROWS][COLS];	        // Buffer for input image
double out_buf[ROWS][COLS];	              // Buffer for output of Laplacian operation, use double to contain negative value
unsigned char out_zero[ROWS][COLS];       // Buffer for zero-crossing result
unsigned char out_Gradient[ROWS][COLS];	  // Buffer for output result of gradiant
unsigned char out_EdgeMap[ROWS][COLS];	  // Buffer for output result of edge map
unsigned char out_hough[ROWS][COLS];      // Buffer for output result of hough transform
unsigned char final_img[ROWS][COLS];      // Buffer for final image with 3 longest lines

// Define the 3x3 Laplacian kernel
int laplacian[3][3] = 
{
  {0, -1, 0},
  {-1,  4, -1},
  {0, -1, 0}
}; 

// This function performs the convolution operation using the Laplacian kernel. 
// It iterates over the image pixels, applies the 3x3 kernel, and stores the result in out_buf. 
// This function also calculates the maximum value of the Laplacian output.
void laplacian_operator() 
{
  // Initialize the max value to 0
  double max = 0;

  // The boundaries are avoided to handle the edges of the image without going out of bounds during convolution.
  // So that starting from 1 and ending at 'ROWS-1'.
  for (int i = 1; i < ROWS - 1; i++) 
  {
    for (int j = 1; j < COLS - 1; j++) 
    {
      // Initializes the sum to store the result of the convolution for each pixel.
      double sum = 0; 

      // Apply the 3x3 Laplacian kernel to the pixel and its neighbors
      for (int k = -1; k <= 1; k++) 
      {
        for (int l = -1; l <= 1; l++)
        {
          sum += in_buf[i + k][j + l] * laplacian[k + 1][l + 1];
        }
      }
      
      // Store the result of the convolution
      out_buf[i][j] = sum; //has detected the spot that matches to the laplacian
    }
  }
}

// This function detects zero-crossings in the Laplacian output and mark them.
// Iterate through each pixel in the image except for the boundaries.
// If a zero-crossing is detected, mark it in 'out_zero'.
void zero_crossing() 
{
  // Initializes a buffer to store the zero-crossing information.
  unsigned char zc_buf[ROWS][COLS] = {0}; 

  // Iterate through each pixel in the image except for the boundaries.
  for (int i = 1; i < ROWS - 1; i++) 
  {
    for (int j = 1; j < COLS - 1; j++) 
    {
      // Store the value of the current pixel in a variable 'center'
      int center = out_buf[i][j]; 
      
      // Check if there is a zero-crossing
      // if the center is positive, and the neighbouring number is negative,
      // mark a zero-crossing
      if ((center > 0 && (out_buf[i + 1][j] < 0 || out_buf[i - 1][j] < 0 || 
                          out_buf[i][j + 1] < 0 || out_buf[i][j - 1] < 0)))
      {
        zc_buf[i][j] = 255; 
      }
    }
  }

  // Copy zero-crossing buffer to output buffer
  for (int i = 0 ; i < ROWS; i++)
  {
    for (int j = 0; j < COLS; j++)
    {
      out_zero[i][j] = zc_buf[i][j];
    }
  }
}

// Applies the Sobel operator on the input image buffer and populates the output image buffers. 
// This function handles the main logic of the edge detection process.
// Calculate Gx and Gy using Sobel masks, and derive the gradient magnitude and edge map.
void sobel_operator() 
{
  // Calculate the boundary for the Sobel operator.
  int M2 = (M - 1) / 2;

  // Define the 3x3 Sobel masks for Gx and Gy.
  double Gx[M][M] = {{-1, 0, 1}, 
                     {-2, 0, 2}, 
                     {-1, 0, 1}};
  
  double Gy[M][M] = {{-1, -2, -1}, 
                     {0, 0, 0},
                     {1, 2, 1}};

  // Iterate through the input image to apply the Sobel operator.
  for (int i = 0; i < ROWS - M2 ; i++) 
  {
    for (int j = 0; j < COLS - M2 ; j++) 
    {
      double sum_Gx = 0;
      double sum_Gy = 0;

      // Convolve the Sobel masks with the image region.
      for (int k = -M2; k <= M2; k++) 
      {
        for (int l = -M2; l <= M2; l++) 
        {
          sum_Gx += in_buf[i + k][j + l] * Gx[k + M2][l + M2];
          sum_Gy += in_buf[i + k][j + l] * Gy[k + M2][l + M2];
        }
      }

      // Calculate temporary values for Gx, Gy, gradient, and edge map.
      //double temp_x = (sum_Gx / (M * M));
      //double temp_y = (sum_Gy / (M * M));
      double temp_gradient = sqrt(sum_Gx * sum_Gx + sum_Gy * sum_Gy);
      // if (out_Gradient > 128 ){out_EdgeMap = 225 }; else{out_EdgeMap = 0}
      // if the gradient at this point is greater than or equal to a pre-specified threshold.
      double temp_edge_map = temp_gradient >= 130 ? 255 : 0;
    
      // Store the results in the respective output buffers.
      //out_Gx[i][j] = (unsigned char)temp_x;
      //out_Gy[i][j] = (unsigned char)temp_y;
      out_Gradient[i][j] = (unsigned char)temp_gradient; // |G(x,y)| = sqrt(Gx^2 + Gy^2)
      //out_EdgeMap[i][j] = temp_edge_map;
      // Generate an edge image where a pixel is an edge point if it is a zero-crossing point and 
      // the gradient at this point is greater than or equal to a pre-specified threshold.
      if (out_zero[i][j] == 255 && temp_edge_map == 255)
      {
        out_EdgeMap[i][j] = 255;
      }
      else
      {
        out_EdgeMap[i][j] = 0;
      }
    }
  }
}

// Apply the Hough transform to the edge map (out_EdgeMap) to detect lines. 
// The accumulator array A is used to store the votes for different rho and theta values. 
// The result is stored in the out_hough buffer.
void hough(unsigned char edge[ROWS][COLS], int A[rMax][tMax], unsigned char hough[ROWS][COLS])
{
  double tRad;  // theta Rad

  // Initialize A[][] to 0's
  for (int r = 0; r < rMax; r++)
  {
    for (int t = 0; t < tMax; t++)
    {
      A[r][t] = 0; 
    }
  }

  // For each edge point (xi,yi) in the edge image f(x,y)
  for (int i = 0; i < ROWS; i++)
  {
    for (int j = 0; j < COLS; j++)
    {
      // (xi,yi) defines a line in θ-ρ parameter space.
      // For (θ = -π/2; θ<=π/2; θ+=Δθ)
      //    ρ = xi * cosθ + yi * sinθ;
      //    A[ρ][θ]++;
      //
      // Check if the current pixel is part of an edge
      if (edge[i][j] == 255)
      {
        // For each possible theta value
        for (int t = 0; t < tMax; t++)
        {
          // Convert theta to radians
          tRad = t * (PI / 180);
          int r = (int)(i * cos(tRad) + j * sin(tRad));

          if (r >= 0 && r < rMax && t >= 0 && t < tMax) 
          {
            // Increment a value in an accumulator array
            A[r][t]++;
          }
        }
      }
    }
  }

  // Find the maximum value in the accumulator array
  int max = 0;
  for (int r = 0; r < rMax; r++)
  {
    for (int t = 0; t < tMax; t++)
    {
      if (A[r][t] > max)
      {
        max = A[r][t];
      }
    }
  }

  // Normalize the accumulator array and convert it to an image
  for (int i = 0; i < ROWS; i++)
  {
    for (int j = 0; j < COLS; j++)
    {
      // Map the rho and theta to i and j
      int rho = (int)(i * (rMax / (float)ROWS));
      int theta = (int)(j * (tMax / (float)COLS));
      // Normalize the value in the range [0, 255]
      hough[i][j] = (unsigned char)(255.0 * A[rho][theta] / max);
    }
  }
}

// Extracts the three longest linear structures from the edge image based on the Hough transform results. 
// The lines are drawn in the out_img buffer.
void drawLines(unsigned char in_img[ROWS][COLS], int A[rMax][tMax], unsigned char out_img[ROWS][COLS])
{
  int n = 3;
  int maxRho[3] = {0, 0, 0};
  int maxTheta[3] = {0, 0, 0};
  int topLine[3] = {0, 0, 0};

  // Find the three most prominent lines in the Hough transform result
  for (int r = 0; r < rMax; r++)
  {
    for (int t = 0; t < tMax; t++)
    {
      int max = A[r][t];

      // Check if the current value is greater than the top three values
      for (int i = 0; i < 3; i++)
      {
        if (max > topLine[i])
        {
          // Shift other values to make room for the new one
          for (int j = 2; j > i; j--)
          {
            topLine[j] = topLine[j - 1];
            maxRho[j] = maxRho[j - 1];
            maxTheta[j] = maxTheta[j - 1];
          }
          // Store the new values
          topLine[i] = max;
          maxRho[i] = r;
          maxTheta[i] = t;
          break;
        }
      }
    }
  }

  // Copy the input image to the output image
  for (int i = 0; i < ROWS; i++)
  {
    for (int j = 0; j < COLS; j++)
    {
      out_img[i][j] = in_img[i][j];
    }
  }

  // Draw the three longest lines on the output image
  for (int i = 0; i < n; i++)
  {
    // Convert Parameters to Image Space
    double rad = maxTheta[i] * PI / 180.0;
    
    for (int x = 0; x < ROWS; x++)
    {
      int y = (maxRho[i] - x * cos(rad)) / sin(rad);
      // Check if the calculated y is within image bounds
      if (y >= 0 && y < COLS)
      {
        out_img[x][y] = 255;
      }
    }
  }
}

int main(int argc, char **argv)
{
  // input file, used to read data from a file and display it
  // output file, used to create files and write data
  FILE *fin, *fout_zero, *fout_gradient, *fout_edge_map, *fout_hough, *fout_final; 
  int n;

  // Check the number of arguments in the command line.
  // If the number of argument is not equal to 3, return error message and exit.
  if (argc != 7) 
  {
    fprintf(stderr, "Usage: %s in.img out_zero.img out_gradient.img out_edge_map.img out_hough.img out_final.img\n", argv[0]);
    exit(1);
  }

  // Open the input image file in the binary mode
  // If cannot open the file, return error message and exit.
  if ((fin = fopen(argv[1], "rb")) == NULL) 
  {
    fprintf(stderr, "ERROR: Cann't open input image file %s\n", argv[1]);
    exit(1);
  }

  // Open the output image file in the binary mode
  // If cannot open the file, return error message and exit.
  if ((fout_zero = fopen(argv[2], "wb")) == NULL) 
  {
    fprintf(stderr, "ERROR: Can't open zero-crossing image file %s\n", argv[2]);
    exit(1);
  }

  if ((fout_gradient = fopen(argv[3], "wb")) == NULL) 
  {
    fprintf(stderr, "ERROR: Can't open gradient image file %s\n", argv[3]);
    exit(1);
  }

  if ((fout_edge_map = fopen(argv[4], "wb")) == NULL) 
  {
    fprintf(stderr, "ERROR: Can't open edge map image file %s\n", argv[4]);
    exit(1);
  }

  if ((fout_hough = fopen(argv[5], "wb")) == NULL) 
  {
    fprintf(stderr, "ERROR: Can't open hough image file %s\n", argv[5]);
    exit(1);
  }

  if ((fout_final = fopen(argv[6], "wb")) == NULL) 
  {
    fprintf(stderr, "ERROR: Can't open final image file %s\n", argv[6]);
    exit(1);
  }

  // Load the input image
  printf("... Load input image\n");
  // read up to ROWS*COLS items of size length from the fin and stores them in the in_buf
  // return the number of full items successfully read
  // If the return is less than count, return error message
  n = fread(in_buf, sizeof(char), ROWS * COLS, fin);
  if (n < ROWS * COLS * sizeof(char))
  {
    fprintf(stderr, "ERROR: Read input image file %s error)\n", argv[1]);
    exit(1);
  }

  // Apply the Laplacian operator and zero-crossing detection
  laplacian_operator();
  zero_crossing();
  // call the sobel_operator function
  sobel_operator();

  int A[rMax][tMax];
  // Apply hough transform to the edge image
  hough(out_EdgeMap, A, out_hough);
  // Extract 3 longest lines
  drawLines(in_buf, A, final_img);

  // Save the output images.
  printf("... Save the output zero-crossing image\n");
  // write up to ROWS*COLS items, each of size bytes in length, from out_buf to fout.
  // If the return is less than count, return error message
  n = fwrite(out_zero, sizeof(char), ROWS * COLS, fout_zero);
  if (n < ROWS * COLS * sizeof(char)) 
  {
    fprintf(stderr, "ERROR: Write output image file %s error)\n", argv[2]);
    exit(1);
  }

  printf("... Save the output gradient image\n");
  n = fwrite(out_Gradient, sizeof(char), ROWS * COLS, fout_gradient);
  if (n < ROWS * COLS * sizeof(char)) 
  {
    fprintf(stderr, "ERROR: Write output image file %s error)\n", argv[3]);
    exit(1);
  }

  printf("... Save the output edge map image\n");
  // write up to ROWS*COLS items, each of size bytes in length, from out_buf to fout.
  // If the return is less than count, return error message
  n = fwrite(out_EdgeMap, sizeof(char), ROWS * COLS, fout_edge_map);
  if (n < ROWS * COLS * sizeof(char)) 
  {
    fprintf(stderr, "ERROR: Write output image file %s error)\n", argv[4]);
    exit(1);
  }

  printf("... Save the output hough image\n");
  // write up to ROWS*COLS items, each of size bytes in length, from out_buf to fout.
  // If the return is less than count, return error message
  n = fwrite(out_hough, sizeof(char), ROWS * COLS, fout_hough);
  if (n < ROWS * COLS * sizeof(char)) 
  {
    fprintf(stderr, "ERROR: Write output image file %s error)\n", argv[5]);
    exit(1);
  }

  printf("... Save the output final image\n");
  // write up to ROWS*COLS items, each of size bytes in length, from out_buf to fout.
  // If the return is less than count, return error message
  n = fwrite(final_img, sizeof(char), ROWS * COLS, fout_final);
  if (n < ROWS * COLS * sizeof(char)) 
  {
    fprintf(stderr, "ERROR: Write output image file %s error)\n", argv[6]);
    exit(1);
  }

  // Close the input file
  fclose(fin);
  // Close the output files
  fclose(fout_zero);
  fclose(fout_gradient);
  fclose(fout_edge_map);
  fclose(fout_hough);
  fclose(fout_final);

  return 0;
}
