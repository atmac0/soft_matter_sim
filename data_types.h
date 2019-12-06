#pragma once

enum Particle_type: int16_t {SQUARE, LONG_RECTANGLE, SHORT_RECTANGLE};

//Define the number of cells, by width and height, of the field.
#define FIELD_WIDTH  1000
#define FIELD_HEIGHT 1000

#define CELL_SIZE 100 //size of each cell in the field in micro-meters. 

#define TIME_INCREMENT 1 //size of the time increment in micro-seconds.

#define SQUARE_SIDE_LENGTH 6000 //please keep this number even. This number is the length of one side of the square in micrometers.
#define NUM_SIDES_SQUARE 4
#define MASS_SQUARE 1 //mass of a square in grams

#define NUM_PARTICLES 2
#define field_t int16_t

struct coord_t
{
  int32_t x;
  int32_t y;
};

struct momentum_t
{
  double x;
  double y;
};

/*********************************/
/* struct png_header		 */
/* {				 */
/*   int length;		 */
/*   int type;			 */
/*   int data;			 */
/*   int CRC;			 */
/* };				 */
/* 				 */
/* const struct png_signature	 */
/* {				 */
/*   char transmission = 137;	 */
/*   char p            = 80;	 */
/*   char n            = 78;	 */
/*   char g            = 71;	 */
/*   char cr           = 13;	 */
/*   char lf           = 10;	 */
/*   char cz           = 26;	 */
/*   char lf2          = 10;     */
/* }				 */
/*********************************/
