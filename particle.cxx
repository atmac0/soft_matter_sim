#include "particle.h"
#include "field.h"

/* PRECONDITION  : None */
/* POSTCONDITION : None */
/* Returns the mass of the particle, based on its particle type. */
uint16_t Particle::get_particle_mass()
{
  if(particle_type == SQUARE)
  {
    return MASS_SQUARE;
  }
}

/* PRECONDITION  : None */
/* POSTCONDITION : A particle has been constructed and all initial values */
/* */
Particle::Particle(Field* field_, field_t par_num, Particle_type par_type, double angular_vel, double orient, double x_vel, double y_vel, coord_t cm_coord)
{
  particle_num       = par_num;
  particle_type      = par_type;

  angular_velocity   = angular_vel;
  orientation        = orient;
  
  x_velocity         = x_vel;
  y_velocity         = y_vel;

  center_mass_coord  = cm_coord;
  field              = field_;

  relative_time      = 0;
  
  if(par_type == SQUARE)
  {
    edge_locations.resize((NUM_SIDES_SQUARE*SQUARE_SIDE_LENGTH/CELL_SIZE)*LINE_THICKNESS); //resize edge location vector to hold all possible edges
    mass = MASS_SQUARE;
    moment_of_inertia = ((double)mass*pow(SQUARE_SIDE_LENGTH,2))/6; //moment of inertia is defined for a flat square as: L = (m*a^2)/6
  }
}  

/* PRECONDITION  : None */
/* POSTCONDITION : The location of an edge has been added the the vector containing the edge locations. */
/* */
void Particle::add_edge_location(coord_t location)
{  
  edge_locations.push_back(location);
}

/* PRECONDITION  : A field has been initialized */
/* POSTCONDITION : The edges of the square have been placed into the field, where the particle number is stored in each cell othat corresponds
                  to an edge point of the square. The location of each edge point placed is stored in a list. Each edge point that is placed 
		  into a cell that already contained an edge point of a different particle will be considered a collision, and stored in a 
		  list that contains all the collision points which will be returned. */
/* The un-rotated orientation of the square is assumed to be as if the top and bottom of the square are parallel with the top and bottom
of the field. A line parallel to these references, pointing to the right, is used as the reference for the x-axis.

........................
.....|-----------|......
.....|...........|......
.....|.....---------->x.
.....|...........|......
.....|-----------|......
........................

The drawing begins at the bottom right corner. This location is found by converting from polar coordinates (using the angle of the orientation),
to cartesian coordinates (using the field as a reference frame), where:

x = rcos(theta)
y = rsin(theta)

where r is calulated by incremented values of an x and y value that are stored relative to the square, using the pythagorean theorem. */
Collisions* Particle::draw_square()
{
  int32_t side_cell_count = SQUARE_SIDE_LENGTH/CELL_SIZE; /*Gives a normalized size of each cell. Aka, if the side length is 20um, and the cell size is 2um,
							the number of cells needed to draw the full edge will be 20um/2um = 10 cells */

  Collisions* collisions_p = new Collisions(SQUARE);
  coord_t coord_field;
  coord_t coord_square;
  for(int32_t j = 0; j < LINE_THICKNESS; j++)
  {
    for(int32_t i = -side_cell_count/2; i <= side_cell_count/2; i++)
    {
      //right side
      coord_square.x = side_cell_count/2 - j;
      coord_square.y = i;
      coord_field = translate_to_field(coord_square);
      field->place_edge_in_field(coord_field, particle_num, collisions_p);
      add_edge_location(coord_field);

      //top side
      coord_square.x = i;
      coord_square.y = -side_cell_count/2 + j;
      coord_field = translate_to_field(coord_square);
      field->place_edge_in_field(coord_field, particle_num, collisions_p);
      add_edge_location(coord_field);
    
      //left side
      coord_square.x = -side_cell_count/2 + j;
      coord_square.y = i;
      coord_field = translate_to_field(coord_square);
      field->place_edge_in_field(coord_field, particle_num, collisions_p);
      add_edge_location(coord_field);
    
      //bottom side
      coord_square.x = i;
      coord_square.y = side_cell_count/2 - j;
      coord_field = translate_to_field(coord_square);
      field->place_edge_in_field(coord_field, particle_num, collisions_p);
      add_edge_location(coord_field);
    }
  }
  return collisions_p;
}

/* PRECONDITION  : None */
/* POSTCONDITION : None */
/* Takes a coordinate of an edge point relative to the center of mass of the particle, and translates it into coordinates relative to the field */
coord_t Particle::translate_to_field(coord_t edge_point)
{
  
  /*From a 2d matrix transform:
    x' = xcos(theta) - ysin(theta)
    y' = xsin(theta) + ycos(theta)
    where x and y are the initial coordinates, and x' and y' are the rotated coordinates */
  int32_t x_prime = edge_point.x*cos(orientation) - edge_point.y*sin(orientation); 
  int32_t y_prime = edge_point.x*sin(orientation) + edge_point.y*cos(orientation);

  // std::cout << "edge_point.x, edge_point.y: (" << edge_point.x <<"," << edge_point.y << ")\n";
  // std::cout << "orientation:                 " << orientation << "\n";
  // std::cout << "x_prime, y_prime: (" << x_prime <<"," << y_prime << ")\n";
  
  //Translate from coordinates relative to the center of mass, to relative to the field
  edge_point.x = center_mass_coord.x + x_prime;
  edge_point.y = center_mass_coord.y + y_prime;	 

  //std::cout << "after: edge_point.x, edge_point.y: (" << edge_point.x <<"," << edge_point.y << ")\n";
  
  //Check if edge is inside of the bounds of the field. If not, have it wrap around.
  if(edge_point.x < 0)
  {
    edge_point.x = FIELD_WIDTH + edge_point.x;
  }
  else
  {
    edge_point.x = edge_point.x%FIELD_WIDTH;
  }

  if(edge_point.y < 0)
  {
    edge_point.y = FIELD_HEIGHT + edge_point.y;
  }
  else
  {
    edge_point.y = edge_point.y%FIELD_HEIGHT;
  }
  
  return edge_point;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* draw edges will identify the particle type, then call the appropriate function to draw the particle */
Collisions* Particle::draw_edges()
{
  if(particle_type == SQUARE)
  {
    return draw_square();
  }
}

/* PRECONDITION  : A field has been initialized */
/* POSTCONDITION : All edge points of the particle be will be cleared and returned to empty cell value. The list containing 
                  the edge locations will also be cleared. */
/* */
void Particle::delete_edges()
{
  for(std::vector<coord_t>::iterator it = edge_locations.begin(); it != edge_locations.end(); it++)
  {
    if(field->get_particle_at(*it) == particle_num)
    {
      field->clear_particle_at(*it);
    }
  }

  edge_locations.clear();  
}

/* PRECONDITION  : None */
/* POSTCONDITION : None */
/* Returns the magnitude of the velocity due to the angular velocity of the furthest edge point from the center of mass of the particle. */
double Particle::find_r_cell_speed()
{
  if(particle_type == SQUARE)
  {
    double r = sqrt(2*pow((SQUARE_SIDE_LENGTH/2)/CELL_SIZE, 2)); //distance from cm to corner
    return abs(r*angular_velocity); //w*r = distance travelled
  }
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
void Particle::propagate()
{
  
  double x_cell_speed = abs(x_velocity)/CELL_SIZE; //convert x velocity to cells per micro-second
  double y_cell_speed = abs(y_velocity)/CELL_SIZE; //convert y velocity to cells per micro-second
  double r_cell_speed = find_r_cell_speed(); //Get rate at which fastest edge tranlates in cells per micro-second
  
  double x_cell_time     = 1/x_cell_speed; //time spent per cell movement
  double y_cell_time     = 1/y_cell_speed;
  double r_cell_time     = 1/r_cell_speed;
  double x_distance_to_travel;
  double y_distance_to_travel;
  
  Collisions* collisions;

  // std::cout << "x_cell_speed: " << x_cell_speed << "\n";
  // std::cout << "y_cell_speed: " << y_cell_speed << "\n";
  // std::cout << "r_cell_speed: " << r_cell_speed << "\n\n";
  
  // std::cout << "x_cell_time: " << x_cell_time << "\n";
  // std::cout << "y_cell_time: " << y_cell_time << "\n";
  // std::cout << "r_cell_time: " << r_cell_time << "\n\n";

  if((r_cell_time > TIME_LIMIT) && (x_cell_time > TIME_LIMIT) && (y_cell_time > TIME_LIMIT))
  {
    rotate_particle(TIME_LIMIT);

    x_distance_to_travel = TIME_LIMIT * abs(x_velocity);
    y_distance_to_travel = TIME_LIMIT * abs(y_velocity);

    translate_y_by_granular(y_distance_to_travel);
    translate_x_by_granular(x_distance_to_travel);

    relative_time += TIME_LIMIT;    
  }
  else if(r_cell_time < x_cell_time && r_cell_time < y_cell_time)
  {
    rotate_particle(r_cell_time);

    x_distance_to_travel = r_cell_time * abs(x_velocity);
    y_distance_to_travel = r_cell_time * abs(y_velocity);
    
    translate_y_by_granular(y_distance_to_travel);
    translate_x_by_granular(x_distance_to_travel);

    relative_time += r_cell_time;
  }
  else if(x_cell_time < y_cell_time)
  {
    rotate_particle(x_cell_time);
    translate_x_by_1();

    y_distance_to_travel = x_cell_time * abs(y_velocity);
    
    translate_y_by_granular(y_distance_to_travel);

    relative_time += x_cell_time;
  }
  else //if(y_cell_time < x_cell_time)
  {
    rotate_particle(y_cell_time);
    translate_y_by_1();

    x_distance_to_travel = y_cell_time * abs(x_velocity);
    
    translate_x_by_granular(x_distance_to_travel);
    
    relative_time += y_cell_time;
  }
}

/* PRECONDITION  : None */
/* POSTCONDITION : The orientation of the particle will be incremented based on the current angular velocity and specified time span. The collisions 
                   incurred due to this rotation will be resolved. */
/* add the rotational displacement (angular velocity * time span) to the orientation. If the orientation goes over 2*PI, take the remainder. */
void Particle::rotate_particle(double time_span)
{
  delete_edges();
  orientation += angular_velocity * time_span;
  orientation = fmod(orientation, M_PI*2);
  
  Collisions* collisions = draw_edges();
  resolve_collisions(collisions);
}


/* PRECONDITION  : None */
/* POSTCONDITION : The previous edges of the particle will be cleared from the field. The position of the center of mass will be incremented
                   or decremented by 1 in the x direction to correspond to the current direction of the x velocity. The new edges corresponding
                   to the new location of the center of mass will be drawn to the field. The collisions incurred from the redrawing will be resolved. */
/* Translate the particle by 1 in the direction appropriate to the current
   velocity of the particle in the x direction. This fuction deletes the current edges, propagates
   the particle in the appropriate direction by 1, then redraws the edges. This
   function returns the collisions that occured from propogation.

   If the particle goes outside of the field, it will wrap around, meaning that if the particle exits
   the bottom, it will reenter thetop, or if it exits the left, it will reenter the right. */
void Particle::translate_x_by_1()
{

  delete_edges();
  
  if(x_velocity > 0)
  {
    center_mass_coord.x++;
    center_mass_coord.x = center_mass_coord.x%FIELD_WIDTH;
  }
  else
  {
    center_mass_coord.x--;
    center_mass_coord.x = center_mass_coord.x%FIELD_WIDTH;
  }

  Collisions* collisions = draw_edges();
  resolve_collisions(collisions);
}

/* PRECONDITION  : None */
/* POSTCONDITION : The previous edges of the particle will be cleared from the field. The position of the center of mass will be incremented
                   or decremented by 1 in the y direction to correspond to the current direction of the y velocity. The new edges corresponding 
		   to the new location of the center of mass will be drawn to the field. The collisions incurred from the redrawing will be resolved. */
/* Translate the particle by 1 in the direction appropriate to the current
   velocity of the particle in the y direction. This fuction deletes the current edges, propagates
   the particle in the appropriate direction by 1, then redraws the edges. This
   function returns the collisions that occured from propogation.
   
   If the particle goes outside of the field, it will wrap around, meaning that if the particle exits
   the bottom, it will reenter thetop, or if it exits the left, it will reenter the right. */
void Particle::translate_y_by_1()
{
  
  delete_edges();
  
  if(y_velocity > 0)
  {
    center_mass_coord.y++;
    center_mass_coord.y = center_mass_coord.y%FIELD_HEIGHT;
  }
  else
  {
    center_mass_coord.y--;
    center_mass_coord.y = center_mass_coord.y%FIELD_HEIGHT;
  }

  Collisions* collisions = draw_edges();
  resolve_collisions(collisions);
}

/* PRECONDITION  : None */
/* POSTCONDITION : The y granular position of the particle will have been incremented or decremeneted by the granularity depending
                   on the direction of the y velocity. If the granular position exceeds the limit of the cell, the center of mass
                   will be translated by 1, and the granular position reset to correspond to the granular position in the new cell.
                   The collisions incurred by a translation will be returned. If the center of mass was not translated, an empty
                   collisions object is returned instead. */
/* */
void Particle::translate_y_by_granular(double granularity)
{
  if(y_velocity > 0)
  {
    y_position_granular += granularity;
  }
  else
  {
    y_position_granular -= granularity;
  }

  if(y_position_granular < 0)
  {
    y_position_granular = 1 + y_position_granular;
    translate_y_by_1();
  }
  if(y_position_granular > 1)
  {
    y_position_granular = y_position_granular - 1;
    translate_y_by_1();
  }
}

/* PRECONDITION  : None */
/* POSTCONDITION : The x granular position of the particle will have been incremented or decremeneted by the granularity depending
                   on the direction of the x velocity. If the granular position exceeds the limit of the cell, the center of mass
                   will be translated by 1, and the granular position reset to correspond to the granular position in the new cell.
                   The collisions incurred by a translation will be returned. If the center of mass was not translated, an empty
                   collisions object is returned instead. */
/* */
void Particle::translate_x_by_granular(double granularity)
{
  if(x_velocity > 0)
  {
    x_position_granular += granularity;
  }
  else
  {
    x_position_granular -= granularity;
  }

  if(x_position_granular < 0)
  {
    x_position_granular = 1 + x_position_granular;
    translate_x_by_1();
  }
  if(x_position_granular > 1)
  {
    x_position_granular = x_position_granular - 1;
    translate_x_by_1();
  }
}

/* PRECONDITION  : None */
/* POSTCONDITION : The linear and angular velocities of the colliding particles will have been altered to reflect the interaction due to
                  all colliding points.*/
/* If the collisions object is empty, the function will return without altering the velocities of the particle. If collisions have occured,
   the particles velocities will be adjusted accordingly. At each collision point, the change in momentum of particle A onto particle B, and
   vice versa, will be calculated seperately. A list storing the total change for each particle will be held until the change for all collisions
   has been calulated. Then, once all the collisions have been resolved, will the change in momentums be converted into change in velocities,
   and the change in velocity will be applied to each corresponding particle.*/
void Particle::resolve_collisions(Collisions* collisions)
{
  
  if(collisions->counts() == 0)
  {
    return;
  }
  
    //dp_t: change in translational momentum; dp_l: change in angular momentum. 1 is the current particle in iteration.
    // 2 is the particle in which the collision happened.
  momentum_t delta_translational[NUM_PARTICLES] = {0};
  double delta_angular[NUM_PARTICLES] = {0};
  momentum_t p1, p2;

  coord_t collision_location;
  field_t particle_1_num, particle_2_num;
  coord_t particle_1_cm_coord, particle_2_cm_coord;

  //Iterate through all the collisions that occured during the most recent drawing
  for(uint32_t i = 0; i < collisions->counts(); i++)
  {   

    //Get the location of the collision and the identifying number of the particles that interacted
    collision_location  = collisions->get_collision_at(i).location;
    particle_1_num      = collisions->get_collision_at(i).particle1;
    particle_2_num      = collisions->get_collision_at(i).particle2;

    //Get the center of mass positions of each particle
    particle_1_cm_coord = particles.at(particle_1_num)->center_mass_coord;
    particle_2_cm_coord = particles.at(particle_2_num)->center_mass_coord;

    //Find the linear momentum of the edge of each particle at the collision location
    p1 = find_linear_momentum_at(collision_location, particle_1_num);
    p2 = find_linear_momentum_at(collision_location, particle_2_num);

    //Find the transfter of momentum due to the individual edge points that collided. The final change will be the sum of the change due to
    //every colliding edge point.
    find_change_in_momentum(particle_1_num, particle_2_num, p1, particle_2_cm_coord, collision_location, delta_translational, delta_angular);
    find_change_in_momentum(particle_2_num, particle_1_num, p2, particle_1_cm_coord, collision_location, delta_translational, delta_angular);
  }

  //This is ~dirty~. I had to use a list of length NUM_PARTICLES because I needed some way to keep track of which particles were interacted with.
  //Most of the elements of the list are going to be empty, however I still iterate through all of them. Implementing some way to keep track of
  //just the particles that were interacted with would be good.
  for(uint32_t i = 0; i < NUM_PARTICLES; i++)
  {
    Particle* particle = particles.at(i);
    //To find the velocity, p=mv, where v=p/m. This is why you must divide by momentum by the mass.
    double average_dx_trans_vel = (delta_translational[i].x/particle->get_mass())/collisions->counts();
    double average_dy_trans_vel = (delta_translational[i].y/particle->get_mass())/collisions->counts();
    double average_dl           = (delta_angular[i]/particle->get_moment_of_inertia())/collisions->counts();
    
    particle->increment_x_velocity(average_dx_trans_vel);
    particle->increment_y_velocity(average_dy_trans_vel);
    particle->increment_angular_velocity(average_dl); 

    // std::cout << "Collision count: " << collisions->counts() << "\n";
    // std::cout << "Particle " << i << ":\n";
    // std::cout << "dx_trans_vel: " << average_dx_trans_vel << "\n";
    // std::cout << "dy_trans_vel: " << average_dy_trans_vel << "\n"; 
    // std::cout << "dl angular  : " << average_dl << "\n\n";       
  }
}

/* PRECONDITION  : d_linear vector points in the direction of a line that intersects with the center of mass. */
/* POSTCONDITION : None */
/*Find if the change in linear momentum is in the direction of the center of mass. This is done by checking if the sign of the
  x or y of the change in linear momentum is the same as the sign as the difference of the coordinate of the center of mass and
  the location of the collision.*/
bool Particle::is_change_in_linear_in_direction_of_cm(coord_t cm_coord, coord_t coll_location, double d_linear_x, double d_linear_y)
{
  int32_t cm_diff_x = cm_coord.x - coll_location.x;
  int32_t cm_diff_y = cm_coord.y - coll_location.y;

  if(((cm_diff_x > 0) && (d_linear_x > 0)) || ((cm_diff_x < 0) && (d_linear_x < 0)))
  {
    return true;
  }
  if(((cm_diff_y > 0) && (d_linear_y > 0)) || ((cm_diff_y < 0) && (d_linear_y < 0)))
  {
    return true;
  }
  
  return false;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/*Find the change in momentum due to particle one onto particle 2. The change in momentum has two parts, the change in translational momentum, and the change in angular momentum. To find the translational change, 
find the projection of the linear momentum vector with the vector pointing from the colliding edge to the center of mass. The x and y components of this
vector will be the */
void Particle::find_change_in_momentum(uint16_t par1, uint16_t par2, momentum_t linear_momentum, coord_t cm_coord, coord_t coll_location, momentum_t dp_t[], double dp_l[])
{
  
  int32_t x_diff      = cm_coord.x - coll_location.x;
  int32_t y_diff      = cm_coord.y - coll_location.y;
  //Correct for edges that have cross the boundary of the field, but not the center of mass.
  if(x_diff < -FIELD_WIDTH/2)
  {
    x_diff -= FIELD_WIDTH;
  }
  if(x_diff > FIELD_WIDTH/2)
  {
    x_diff += FIELD_WIDTH;
  }
  if(y_diff < -FIELD_HEIGHT/2)
  {
    y_diff -= FIELD_HEIGHT;
  }
  if(y_diff > FIELD_HEIGHT/2)
  {
    y_diff += FIELD_HEIGHT;
  }
  
  double   r_mag       = sqrt(pow(x_diff,2) + pow(y_diff,2));
  double   lin_mom_mag = sqrt(pow(linear_momentum.x,2) + pow(linear_momentum.y,2));  

  double d_linear_x = (linear_momentum.x * abs(x_diff))/r_mag;
  double d_linear_y = (linear_momentum.y * abs(y_diff))/r_mag;

  if(is_change_in_linear_in_direction_of_cm(cm_coord, coll_location, d_linear_x, d_linear_y))
  {
    dp_t[par2].x += d_linear_x;
    dp_t[par2].y += d_linear_y;
    
    dp_t[par1].x -= d_linear_x;
    dp_t[par1].y -= d_linear_y;
  }
  
  double dp_lx = (linear_momentum.x *  abs(y_diff))/r_mag;
  double dp_ly = (linear_momentum.y * -abs(x_diff))/r_mag;

  int32_t l_direction = determine_direction_of_angular_change(x_diff, y_diff, dp_lx, dp_ly);

  double d_angular_momentum = l_direction*sqrt( pow(dp_lx,2) + pow(dp_ly,2) );
  
  dp_l[par2] += d_angular_momentum;
  dp_l[par1] -= d_angular_momentum;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/*Returns -1 if the angular change is in the clockwise direction, 1 if in the counter clockwise direction. Set up a matrix where the elements are
organzied as so:

x_diff, y_diff
dp_lx , dp_ly

if the direction of change is in the couter clockwise direction, cross multiplications shows: x_diff*dp_ly >= 0 && dp_lx * y_diff <= 0 */
int32_t Particle::determine_direction_of_angular_change(int32_t x_diff, int32_t y_diff, double dp_lx, double dp_ly)
{
  if(((x_diff * dp_ly) >= 0) && ((dp_lx * y_diff) <= 0))
  {
    return -1;
  }

  return 1;
  
  /*for counter clockwise:
    diff: +x, 0y
    p   : 0x, +y

    diff: +x, +y
    p:    -x, +y

    diff: 0x, +y
    p   : -x, 0y

    diff: -x, +y
    p   : -x, -y

    diff: -x, 0y
    p   : 0x, -y

    diff: -x, -y
    p   : +x, -y

    diff: 0x, -y
    p   : +x, 0y

    diff: +x, -y
    p   : +x, +y    
  */
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/*Find the linear momentum at a point on the edge of a particle.
Notation:
L: angular momentum
v: translational velocity
r: disance from point of rotation (center of mass)
w: angular velocity
m: mass of particle
I: moment of inertia

L = rp; L=Iw; => p = Iw/r

These relations are used to find the magnitude of the momentum at a point on the edge. The components
are then found using */
momentum_t Particle::find_linear_momentum_at(coord_t point, field_t particle_num)
{

  Particle* particle = particles.at(particle_num);
  momentum_t linear_momentum;

  //Find the components due to translational momentum
  linear_momentum.x = particle->get_x_velocity()*particle->get_mass();
  linear_momentum.y = particle->get_y_velocity()*particle->get_mass();

  //Find the components due to angular momentum
  /*  To find the linear momentum due to the angular velocity, imagine you draw a circle about the center of mass where
  the edge intersects with the point at which a collision happened. Draw a line from the center of mass to this point.
  The direction of the momentum is the negative inverse of the slope of this line. However, the direction is opposite if 
  the particle is rotating counter clockwise.

  ................+.......
  .....|-----------*+..... The stars represent the r vector. The + represent the perpendicular momentum vector.
  .....|........*..|..+....
  .....|.....*--------->x.
  .....|...........|......
  .....|-----------|......
  ........................
*/

  uint32_t x_diff = point.x - particle->get_center_mass_coord().x;
  uint32_t y_diff = point.y - particle->get_center_mass_coord().y;

  if(x_diff < -FIELD_WIDTH/2)
  {
    x_diff -= FIELD_WIDTH;
  }
  if(x_diff > FIELD_WIDTH/2)
  {
    x_diff += FIELD_WIDTH;
  }
  if(y_diff < -FIELD_HEIGHT/2)
  {
    y_diff -= FIELD_HEIGHT;
  }
  if(y_diff > FIELD_HEIGHT/2)
  {
    y_diff += FIELD_HEIGHT;
  }
  
  double r        = sqrt(pow(x_diff,2) + pow(y_diff,2));
  
  double w        = particle->get_angular_velocity();
  double I        = particle->get_moment_of_inertia();
  //The linear momentum due to the angular momentum at a distance r from the square
  double p_from_ang = I*w/r;

  //Check if the difference in the x or y is zero. If it is, all of the angular momentum is linear in the corresponding x or y direction.
  if(x_diff == 0)
  {
    linear_momentum.x += p_from_ang;
    return linear_momentum;
  }
  if(y_diff == 0)
  {
    linear_momentum.y += p_from_ang;
    return linear_momentum;
  }

  //slope of the line from the center of mass to the point
  double slope = y_diff/x_diff;
  //The slope of the line in the direction of the angular momentum
  double p_ang_slope = 1/slope;   
  double theta = atan(p_ang_slope);
  
  linear_momentum.x += p_from_ang*sin(theta);
  linear_momentum.y += p_from_ang*cos(theta);
  return linear_momentum;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
uint32_t Particle::get_mass()
{
  return mass;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
void Particle::increment_y_velocity(double dy)
{
  y_velocity += dy;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
void Particle::increment_x_velocity(double dx)
{
  x_velocity += dx;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
void Particle::increment_angular_velocity(double dw)
{
  angular_velocity += dw;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
double Particle::get_y_velocity()
{
  return y_velocity;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
double Particle::get_x_velocity()
{
  return x_velocity;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
coord_t Particle::get_center_mass_coord()
{
  return center_mass_coord;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
double Particle::get_angular_velocity()
{
  return angular_velocity;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
void Particle::set_particles_vector(std::vector<Particle*>* par_vec)
{
  particles = *par_vec;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
double Particle::get_moment_of_inertia()
{
  return moment_of_inertia;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
double Particle::get_relative_time()
{
  return relative_time;
}

/* PRECONDITION  : */
/* POSTCONDITION : */
/* */
field_t Particle::get_particle_num()
{
  return particle_num;
}
