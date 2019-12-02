#pragma once

enum Particle_type: int16_t {SQUARE, LONG_RECTANGLE, SHORT_RECTANGLE};

//Define the number of cells, by width and height, of the field.
#define FIELD_WIDTH  100
#define FIELD_HEIGHT 100

#define CELL_SIZE 1 //size of each cell in the field in micro-meters. 

#define TIME_INCREMENT 1 //size of the time increment in micro-seconds.

#define SQUARE_SIDE_LENGTH 20 //please keep this number even. This number is the length of one side of the square in micrometers.
#define NUM_SIDES_SQUARE 4

struct coord_t
{
  int32_t x;
  int32_t y;
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
