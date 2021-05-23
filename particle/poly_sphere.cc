#include "poly_sphere.h"

double part_poly_sphere::r_coordinate_box(int j) { return coordinate_box[j]; }

vector <double> part_poly_sphere::r_coordinate_box_vec() { return coordinate_box; }

double part_poly_sphere::r_coordinate_true(int j) { return coordinate_true[j]; }

vector <double> part_poly_sphere::r_coordinate_true_vec() { return coordinate_true; }

double part_poly_sphere::r_coordinate_affine(int j) { return coordinate_affine[j]; }

vector <double> part_poly_sphere::r_coordinate_affine_vec() { return coordinate_affine; }

int part_poly_sphere::r_type() { return type; }

double part_poly_sphere::r_size() { return size; }

void part_poly_sphere::c_coordinate_box(int j, double new_coordinate_box) { coordinate_box[j] = new_coordinate_box; }

void part_poly_sphere::c_coordinate_true(int j, double new_coordinate_true) { coordinate_true[j] = new_coordinate_true; }

void part_poly_sphere::c_coordinate_affine(int j, double new_coordinate_affine) { coordinate_affine[j] = new_coordinate_affine; }

void part_poly_sphere::c_type(int new_type) { type = new_type; }

void part_poly_sphere::c_size(double new_size) { size = new_size; }
