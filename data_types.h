#pragma once

#define field_t int16_t

enum Particle_type        : int16_t {SQUARE, LONG_RECTANGLE, SHORT_RECTANGLE};
enum Lowest_particle_time : field_t {X_TIME, Y_TIME, ROTATIONAL_TIME};
enum Square_side          : int16_t {RIGHT, TOP, LEFT, BOTTOM};

//Define the number of cells, by width and height, of the field.
#define FIELD_WIDTH  1000
#define FIELD_HEIGHT 1000

#define CELL_SIZE 1 //size of each cell in the field in micro-meters. 

#define TIME_INCREMENT 1 //size of the time increment in micro-seconds.
const double TIME_LIMIT = 0.1; //maximum length of time 

#define LINE_THICKNESS 4

#define SQUARE_SIDE_LENGTH 20 //please keep this number even. This number is the length of one side of the square in micrometers.
#define NUM_SIDES_SQUARE 4
#define MASS_SQUARE 1 //mass of a square in grams

<<<<<<< HEAD
#define LINE_THICKNESS 6

#define NUM_PARTICLES 20

=======
#define NUM_PARTICLES 3
#define field_t int16_t
>>>>>>> parent of 7683a69... Fixed collision resolution bug where particles dragged eachother along. Added documentation to particle class. Cleaned up center of mass translation functions.

struct coord_t
{
  field_t x;
  field_t y;
};

struct momentum_t
{
  double x;
  double y;
};
