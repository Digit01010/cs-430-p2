#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int line = 1;

// Holds an rgb triple of a pixel
typedef struct Pixel {
  unsigned char red, green, blue;
} Pixel;

// Holds information about the header of a ppm file
typedef struct Header {
   unsigned char magicNumber;
   unsigned int width, height, maxColor;
} Header;

// Plymorphism in C

typedef struct {
  int kind; // 0 = camera, 1 = sphere, 2 = plane
  double color[3];
  union {
    struct {
      double width;
      double height;
    } camera;
    struct {
      double position[3];
      double radius;
    } sphere;
    struct {
      double position[3];
      double normal[3];
    } plane;
  };
} Object;

double sphere_intersection(double*, double*, double*, double);
double plane_intersection(double*, double*, double*, double*);
void writeP3(Pixel *, Header, FILE *);
int next_c(FILE*);
void expect_c(FILE*, int);
void skip_ws(FILE*);
char* next_string(FILE*);
double next_number(FILE*);
double* next_vector(FILE*);
Object** read_scene(char*);

static inline double sqr(double v) {
  return v*v;
}


static inline void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "Error: Incorrect number of arguments.\n");
    printf("Usage: raycast width height input.json output.ppm\n");
    return(1);
  }
  
  int N = atoi(argv[1]);
  int M = atoi(argv[2]);

  Object** objects = read_scene(argv[3]);
  /*Object** objects;
  objects = malloc(sizeof(Object*)*3);
  objects[0] = malloc(sizeof(Object));
  objects[0]->kind = 1;
  objects[0]->color[0] = 1;
  objects[0]->color[1] = 0;
  objects[0]->color[2] = 0;
  objects[0]->sphere.radius = 2;
  // object[0]->teapot.handle_length = 2;
  objects[0]->sphere.position[0] = 1;
  objects[0]->sphere.position[1] = 1;
  objects[0]->sphere.position[2] = 10;
  
  objects[1] = malloc(sizeof(Object));
  objects[1]->kind = 2;
  objects[1]->color[0] = 0;
  objects[1]->color[1] = 1;
  objects[1]->color[2] = 0;
  // object[0]->teapot.handle_length = 2;
  objects[1]->plane.position[0] = 0;
  objects[1]->plane.position[1] = -1;
  objects[1]->plane.position[2] = 0;
  objects[1]->plane.normal[0] = 0;
  objects[1]->plane.normal[1] = 1;
  objects[1]->plane.normal[2] = 0;
  
  objects[2] = NULL;*/
  
  double cx = 0;
  double cy = 0;
  double h = 1;
  double w = 1;
  
  Pixel *buffer = malloc(sizeof(Pixel) * N * M);
  
  double pixheight = h / M;
  double pixwidth = w / N;
  for (int y = M; y > 0; y -= 1) {
    for (int x = 0; x < N; x += 1) {
      double Ro[3] = {0, 0, 0};
      // Rd = normalize(P - Ro)
      double Rd[3] = {
        cx - (w/2) + pixwidth * (x + 0.5),
        cy - (h/2) + pixheight * (y + 0.5),
        1
      };
      normalize(Rd);

      double best_t = INFINITY;
      int best_i = 0;
      for (int i=0; objects[i] != 0; i += 1) {
        double t = 0;

        switch(objects[i]->kind) {
        case 0:
          t = -1;
          break;
        case 1:
          t = sphere_intersection(Ro, Rd,
                                    objects[i]->sphere.position,
                                    objects[i]->sphere.radius);
          break;
        case 2:
          t = plane_intersection(Ro, Rd,
                                    objects[i]->plane.position,
                                    objects[i]->plane.normal);
          break;
        default:
          fprintf(stderr, "Error: Programmer forgot to implement an intersection.");
          exit(1);
        }
        if (t > 0 && t < best_t) {
          best_t = t;
          best_i = i;
        }
      }
      int p = (M - y)*N + x;
      if (best_t > 0 && best_t != INFINITY) {
        buffer[p].red = (char) objects[best_i]->color[0] * 255;
        buffer[p].green = (char) objects[best_i]->color[1] * 255;
        buffer[p].blue = (char) objects[best_i]->color[2] * 255;
      } else {
        buffer[p].red = 0;
        buffer[p].green = 0;
        buffer[p].blue = 0;
      }
      
    }
  }
  
  FILE* output = fopen(argv[4], "w");
  
  if (output == NULL) {
    fprintf(stderr, "Error: Could not write to file \"%s\"\n", filename);
    exit(1);
  }
  
  Header outHeader;
  outHeader.magicNumber = 3;
  outHeader.maxColor = 255;
  outHeader.width = N;
  outHeader.height = M;
  
  writeP3(buffer, outHeader, output);
  
  
  return 0;
}

double sphere_intersection(double* Ro, double* Rd,
                             double* C, double r) {
  double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
  double b = 2 * (Rd[0] * (Ro[0] - C[0]) + Rd[1] * (Ro[1] - C[1]) + Rd[2] * (Ro[2] - C[2]));  
  double c = sqr(Ro[0] - C[0]) + sqr(Ro[1] - C[1]) + sqr(Ro[2] - C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);
  
  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}

double plane_intersection(double* Ro, double* Rd,
                             double* P, double* n) {
  double t = -(n[0]*(Ro[0]-P[0]) + n[1]*(Ro[1]-P[1]) + n[2]*(Ro[2]-P[2])) / (n[0]*Rd[0] + n[1]*Rd[1] + n[2]*Rd[2]);
  return t;
}


/*double cylinder_intersection(double* Ro, double* Rd,
                             double* C, double r) {
  // Step 1. Find the equation for the object you are
  // interested in..  (e.g., cylinder)
  //
  // x^2 + z^2 = r^2
  //
  // Step 2. Parameterize the equation with a center point
  // if needed
  //
  // (x-Cx)^2 + (z-Cz)^2 = r^2
  //
  // Step 3. Substitute the eq for a ray into our object
  // equation.
  //
  // (Rox + t*Rdx - Cx)^2 + (Roz + t*Rdz - Cz)^2 - r^2 = 0
  //
  // Step 4. Solve for t.
  //
  // Step 4a. Rewrite the equation (flatten).
  //
  // -r^2 +
  // t^2 * Rdx^2 +
  // t^2 * Rdz^2 +
  // 2*t * Rox * Rdx -
  // 2*t * Rdx * Cx +
  // 2*t * Roz * Rdz -
  // 2*t * Rdz * Cz +
  // Rox^2 -
  // 2*Rox*Cx +
  // Cx^2 +
  // Roz^2 -
  // 2*Roz*Cz +
  // Cz^2 = 0
  //
  // Steb 4b. Rewrite the equation in terms of t.
  //
  // t^2 * (Rdx^2 + Rdz^2) +
  // t * (2 * (Rox * Rdx - Rdx * Cx + Roz * Rdz - Rdz * Cz)) +
  // Rox^2 - 2*Rox*Cx + Cx^2 + Roz^2 - 2*Roz*Cz + Cz^2 - r^2 = 0
  //
  // Use the quadratic equation to solve for t..
  double a = (sqr(Rd[0]) + sqr(Rd[2]));
  double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[2] * Rd[2] - Rd[2] * C[2]));
  double c = sqr(Ro[0]) - 2*Ro[0]*C[0] + sqr(C[0]) + sqr(Ro[2]) - 2*Ro[2]*C[2] + sqr(C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);
  
  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}
*/

// Writes P3 data
void writeP3(Pixel *buffer, Header h, FILE *fh) {
  // Write the header
  fprintf(fh, "P%d\n%d %d\n%d\n", h.magicNumber, h.width, h.height, h.maxColor);
  // Write the ascii data
  for (int i = 0; i < h.width * h.height; i++) {
     fprintf(fh, "%d\n%d\n%d\n", buffer[i].red, buffer[i].green, buffer[i].blue);
  }
}

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);    
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }  
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);      
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);      
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  if (ferror(json)) {
    fprintf(stderr, "Error: Could not read number on line %d.\n", line);
    exit(1);
  }
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}


Object** read_scene(char* filename) {
  int c;
  
  Object** objects;
  objects = malloc(sizeof(Object*)*129);
  
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }
  
  skip_ws(json);
  
  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects
  int objcnt = 0;
  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return NULL;
    }
    if (c == '{') {
      skip_ws(json);
    
      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
        fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
        exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);
      if (strcmp(value, "camera") == 0) {     
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 0;
      } else if (strcmp(value, "sphere") == 0) {
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 1;
      } else if (strcmp(value, "plane") == 0) {
        objects[objcnt] = malloc(sizeof(Object));
        objects[objcnt]->kind = 2;
      } else {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        exit(1);
      }

      skip_ws(json);

      while (1) {
        // , }
        c = next_c(json);
        if (c == '}') {
          // stop parsing this object
          objcnt += 1;
          if (objcnt > 128) {
            exit(1);
          }
          break;
        } else if (c == ',') {
          // read another field
          skip_ws(json);
          char* key = next_string(json);
          skip_ws(json);
          expect_c(json, ':');
          skip_ws(json);
          
          if (strcmp(key, "width") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 0:
              objects[objcnt]->camera.width = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "height") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 0:
              objects[objcnt]->camera.height = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "radius") == 0) {
            double value = next_number(json);
            switch (objects[objcnt]->kind) {
            case 1:
              objects[objcnt]->sphere.radius = value;
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "color") == 0) {
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 0:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            default:
              objects[objcnt]->color[0] = value[0];
              objects[objcnt]->color[1] = value[1];
              objects[objcnt]->color[2] = value[2];
              break;
            }
          } else if (strcmp(key, "position") == 0){
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 1:
              objects[objcnt]->sphere.position[0] = value[0];
              objects[objcnt]->sphere.position[1] = value[1];
              objects[objcnt]->sphere.position[2] = value[2];
              break;
            case 2:
              objects[objcnt]->plane.position[0] = value[0];
              objects[objcnt]->plane.position[1] = value[1];
              objects[objcnt]->plane.position[2] = value[2];
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else if (strcmp(key, "normal") == 0) {
            double* value = next_vector(json);
            switch (objects[objcnt]->kind) {
            case 2:
              objects[objcnt]->plane.normal[0] = value[0];
              objects[objcnt]->plane.normal[1] = value[1];
              objects[objcnt]->plane.normal[2] = value[2];
              break;
            default:
              fprintf(stderr, "Error: Unexpected key on line %d.\n", line);
              exit(1);
              break;
            }
          } else {
            fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
                    key, line);
            //char* value = next_string(json);
          }
          skip_ws(json);
        } else {
          fprintf(stderr, "Error: Unexpected value on line %d\n", line);
          exit(1);
        }
      }
      skip_ws(json);
      c = next_c(json);
      if (c == ',') {
        // noop
        skip_ws(json);
      } else if (c == ']') {
        objects[objcnt] = NULL;
        fclose(json);
        return objects;
      } else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    }

  }

  
  return NULL;
}

