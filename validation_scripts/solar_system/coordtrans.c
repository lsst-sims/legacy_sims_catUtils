#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define skipline(ptr) while((getc(ptr)!='\n') && !(feof(ptr)))

void print_help(void);
/* this is the one that is used to rotate coords */
void radectransform(float inc, float ra, float dec,
		     float *lat, float *lon);
/* i don't actually use the radecfancy anymore, it was an experiment */
void radecfancy(float node, float inc, float ra, float dec,
		float *lat, float *lon) ;


void print_help() {
  printf("transform [-b filename] [-r, e, i, k]  RA DEC\n");
  printf("  tranforms a given ra and dec into ecliptic plane\n");
  printf("  coordinates and invariable plane coordinates\n");
  printf("  or vice versa\n");
  printf("[-b file]: input from file of ra/dec - must include filename\n");
  printf("[-r,e,i,k] : use to indicate type of incoming coordinates\n");
  printf("   note: 'k' coordinates are KBO plane coords, as by MB\n");
  printf("Coordinates should be given in 'normal' units, J2000\n");
  printf("and long(or ra) then lat (or dec)\n");
  exit(-1);
}

float deg_per_hrs = 15;
float rad_per_deg = 0.0174533;
float deg_per_rad = 57.2957795130823;

int main(int argc, char *argv[]) 
{
  int rah, ramin;
  float ras, ra;
  int  decd, decmin;
  float decs, dec;
  float eclon, eclat;
  float invlon, invlat;
  float kbolon, kbolat;
  float ecinc, invinc, ecnode, invnode;
  float invecnode, invecinc;
  float kbonode, kboinc;
  int i = 1;
  int sign = 1;
  float temp;
  char stemp[1];
  FILE *fptr;
  
  /* just do everything in J2000 and we should be fine */
  ecnode = 0; /* relative ra/dec */
  ecinc = (23.439291)*rad_per_deg;
  invnode = (3+51.0/60.0+9.4/3600.0)*rad_per_deg; /* relative ra/dec */
  invinc = (23+32.0/3600.0)*rad_per_deg;
  invecnode = (107+34.0/60.0+57.7/3600)*rad_per_deg; /*relative to ec */
  invecinc = (1+34.0/60.0+43.3/3600)*rad_per_deg;
  kbonode = 81.6 * rad_per_deg; /* relative to ecliptic */
  kboinc = 1.86 * rad_per_deg; /* relative to ecliptic */
  /* relative to J2000 coordinates! */
  /* printf("ecliptic node and inclination: %f %f\n", ecnode, ecinc); */
  /* printf("invariable node and inclination: %f %f\n", invnode, invinc); */
  
  if (argc < 3) 
    print_help();

  /* read in options */
  if (argv[i][0] == '-' && atof(argv[i])==0) { /*specifying an option */
    /* want to transform ra/dec to others */
    if (strcasecmp(argv[i]+1, "r")==0) { /* giving ra/dec */
      i++;
      /* read ra and dec */
      if (argc == 4) { /* radec must be given in hrs/deg only */
	ra = atof(argv[i++]);
	dec = atof(argv[i++]);
      }
      /* you got two options - simple hrs/deg only or the full thing */
      else {
	rah = atoi(argv[i++]);
	ramin = atoi(argv[i++]);
	ras = atof(argv[i++]);
	ra = rah + ramin/60.0 + ras/3600.0;
	decd = atoi(argv[i++]);
	decmin = atoi(argv[i++]);
	decs = atof(argv[i++]);
	if (decd < 0) { /*dec less than zero */
	  sign = -1;
	  decd = -decd;
	}
	/* need to watch for this special case of dec = -0 x x */
	if (decd == 0 && argv[i-3][0] == '-') sign = -1;
	dec = sign*(decd + decmin/60.0 + decs/3600.0);
      }
      printf("RA and Dec are %g (hrs) %g (deg) \n", ra, dec);
  
      ra = ra*deg_per_hrs;
      ra = ra*rad_per_deg;
      dec = dec*rad_per_deg;

      printf("RA and Dec in radians %g %g\n", ra,dec);

      /* transform radec into ecliptic coords */
      radectransform(ecinc, ra, dec, &eclon, &eclat);
      printf("Ecliptic lon and lat %g %g\n", eclon, eclat);
      /* transform radec into invariable plane coords- rel to radec */
      radectransform(invinc, (ra-invnode), dec, &invlon, &invlat);
      printf("Invariable lon and lat %g %g\n", invlon, invlat);
      /* transform radec into kboplane coords - relative to ecliptic */
      eclon = eclon*rad_per_deg;
      eclat = eclat*rad_per_deg;
      radectransform(kboinc, (eclon-kbonode), eclat, &kbolon, &kbolat);
      printf("KBOplane lon and lat %g %g\n", kbolon, kbolat);

      /*radecfancy(invnode, invinc, ra, dec, &invlon, &invlat);  */
      /*printf("Invariable lon and lat %g %g\n", invlon, invlat);  */
      

    } /* done with everything for ra/dec input */

    /* or maybe you are inputting ecliptic coordinates*/
    else if (strcasecmp(argv[i]+1, "e")==0) {
      /* ecliptic coordinate input */
      i++;
      /* read ec lat/long in degrees */
      eclon = atof(argv[i++]);
      eclat = atof(argv[i++]);

      eclat = eclat*rad_per_deg;
      eclon = eclon*rad_per_deg;
      printf("Ecliptic lon/lat : %f %f\n", eclon, eclat);

      /* transform eclonlat into RADEC */
      radectransform(-ecinc, eclon, eclat, &ra, &dec);
      ra = ra / deg_per_hrs;
      printf("RA (hrs) and dec (deg) : %f %f\n", ra, dec);
      rah = (int) ra;
      ramin = (int) ((ra-rah)*60);
      ras = ((ra-rah)*60 - ramin)*60;
      if (dec>0) {
	decd = (int) dec;
	decmin = (int) ((dec-decd)*60);
	decs = ((dec-decd)*60 - decmin)*60;
	printf("RA: %d %d %.2f DEC: %d %d %.2f\n", 
	       rah, ramin, ras, decd, decmin, decs);
      }
      else {
	dec = -dec;
	decd = (int) dec;
	decmin = (int) ((dec-decd)*60);
	decs = ((dec-decd)*60 - decmin)*60;
	printf("RA: %d %d %.2f DEC: -%d %d %.2f\n", 
	       rah, ramin, ras, decd, decmin, decs);
      }
      /* transform eclonlat into invariable coords, wrt ecliptic plane */
      radectransform(invecinc, (eclon-invecnode), eclat, &invlon, &invlat);
      printf("Invariable lon/lat wrt ec : %f %f \n", invlon, invlat);
      ra = ra*deg_per_hrs*rad_per_deg;
      dec = dec*rad_per_deg;
      radectransform(invinc, (ra-invnode), dec, &invlon, &invlat);
      printf("Invariable lon/lat wrt radec : %f %f \n", invlon, invlat);
      /* transform eclonlat into kbolonlat */
      radectransform(kboinc, (eclon-kbonode), eclat, &kbolon, &kbolat);
      printf("KBOplane lon/lat wrt ecliptic %f %f\n", kbolon, kbolat);
      
    }
    
    else if (strcasecmp(argv[i]+1,"i") ==0) { /* invariable plane input */
      i++;
      invlon = atof(argv[i++]);
      invlat = atof(argv[i++]);
      invlat = invlat*rad_per_deg;
      invlon = invlon*rad_per_deg;
      printf("Invariable lon/lat, WRT RA/DEC (rad) : %f %f\n", invlon, invlat);

      /* transform inv lonlat into radec */
      radectransform(-invinc, invlon, invlat, &ra, &dec);
      ra = ra+invnode*deg_per_rad;
      ra = ra / deg_per_hrs;
      if (ra>24) ra=ra-24;
      printf("RA (hrs) and dec (deg) : %f %f\n", ra, dec);
      rah = (int) ra;
      ramin = (int) ((ra-rah)*60);
      ras = ((ra-rah)*60 - ramin)*60;
      if (dec>0) {
	decd = (int) dec;
	decmin = (int) ((dec-decd)*60);
	decs = ((dec-decd)*60 - decmin)*60;
	printf("RA: %d %d %.2f DEC: %d %d %.2f\n", 
	       rah, ramin, ras, decd, decmin, decs);
      }
      else {
	dec = -dec;
	decd = (int) dec;
	decmin = (int) ((dec-decd)*60);
	decs = ((dec-decd)*60 - decmin)*60;
	printf("RA: %d %d %.2f DEC: -%d %d %.2f\n", 
	       rah, ramin, ras, decd, decmin, decs);
      }
      /* transform in inv lon/lat into ecliptic */
      radectransform(-invecinc, invlon, invlat, &eclon, &eclat);
      eclon = eclon+invecnode*deg_per_rad;
      if (eclon > 360) eclon=eclon-360;
      printf("Ecliptic lon/lat wrt invariable : %f %f \n", eclon, eclat);
      ra = ra*deg_per_hrs*rad_per_deg;
      dec = dec*rad_per_deg;
      /* transform inv lon/lat into ecliptic, using inv wrt radec */
      radectransform(ecinc, ra, dec, &eclon, &eclat);
      printf("Ecliptic lon/lat wrt radec : %f %f \n", eclon, eclat);
      /* for (temp = 0; temp<7; temp=temp+0.02) { */
      /*  for (invlat = -0.2; invlat < 0.3; invlat = invlat+0.2) { */
      /*  radectransform(-invecinc, temp, invlat*rad_per_deg, &eclon, &eclat);*/
      /* eclon = eclon + invecnode*deg_per_rad; */
      /*if (eclon>360) eclon=eclon-360; */
      /*printf("%f %f %f\n", temp, eclon, eclat); */
      /* transform inv lon/lat into KBO plane, using ecliptic coods */
      eclon = eclon*rad_per_deg;
      eclat = eclat*rad_per_deg;
      radectransform(kboinc, (eclon-kbonode), eclat, &kbolon, &kbolat);
      printf("KBOplane lon/lat %f %f (deg)\n", kbolon, kbolat);
    }
    
    else if (strcasecmp(argv[i]+1, "k")==0) { /* transform from kbo plane */
      i++;
      kbolon = atof(argv[i++]);
      kbolat = atof(argv[i++]);
      kbolat = kbolat*rad_per_deg;
      kbolon = kbolon*rad_per_deg;
      printf("KBOplane lon/lat is %f %f (rad)\n", kbolon, kbolat);
      /* transform kbolon/lat into ecliptic coords */
      radectransform(-kboinc, kbolon, kbolat, &eclon, &eclat);
      eclon = eclon + kbonode*deg_per_rad;
      if (eclon>360) eclon = eclon-360.;
      printf("Ecliptic latitude/long (deg) %f %f\n", eclon, eclat);
      /* transform kbolon/lat into radec .. use ecliptic */
      radectransform(-ecinc, eclon, eclat, &ra, &dec);
      ra = ra / deg_per_hrs;
      if (ra>24) ra=ra-24;
      printf("RA (hrs) and dec (deg) : %f %f\n", ra, dec);
      rah = (int) ra;
      ramin = (int) ((ra-rah)*60);
      ras = ((ra-rah)*60 - ramin)*60;
      if (dec>0) {
	decd = (int) dec;
	decmin = (int) ((dec-decd)*60);
	decs = ((dec-decd)*60 - decmin)*60;
	printf("RA: %d %d %.2f DEC: %d %d %.2f\n", 
	       rah, ramin, ras, decd, decmin, decs);
      }
      else {
	dec = -dec;
	decd = (int) dec;
	decmin = (int) ((dec-decd)*60);
	decs = ((dec-decd)*60 - decmin)*60;
	printf("RA: %d %d %.2f DEC: -%d %d %.2f\n", 
	       rah, ramin, ras, decd, decmin, decs);
      }
      /* for (temp = 0; temp<7; temp=temp+0.02) { */
      /*  for (kbolat = -0.2; kbolat < 0.3; kbolat = kbolat+0.2) { */
      /*  radectransform(-invecinc, temp, invlat*rad_per_deg, &eclon, &eclat);*/
      /* eclon = eclon + invecnode*deg_per_rad; */
      /*if (eclon>360) eclon=eclon-360; */
      /*printf("%f %f %f\n", temp, eclon, eclat); */
    }



    /* batch mode input - from file */
    else if (strcasecmp(argv[i]+1, "b")==0) { 
      /* input from radec file */
      i++;
      if ((fptr = fopen(argv[i], "r"))==NULL) {
	printf("Can't open file\n");
	print_help();
      }
      fprintf(stderr, "Please be sure that any declinations of the form -00 xx xx are edited to be 00 -xx xx\n");
      fprintf(stderr, "  -- sorry, that's just the way it has to be\n");
      printf("#RA\tDec\tEclon\tEclat\tInvlon\tInvlat\n");
      while (!feof(fptr)) {
	if (fscanf(fptr, "%d %d %f %d %d %f", 
		   &rah, &ramin, &ras, &decd, &decmin, &decs) <6)
	  break;
	skipline(fptr);
	ra = rah + ramin/60.0 + ras/3600.0;
	ra = ra*deg_per_hrs*rad_per_deg;
	sign=1;
	if (decd<0) { 
	  sign = -1;
	  decd = -decd;
	}
	/* need to watch out for special case of -00 xx xx for dec */
	/* will have to edit original file for kludge */
	if (decmin<0) {
	  sign=-1;
	  decmin = -decmin;
	}
	dec = sign*(decd + decmin/ 60.0 + decs/3600.0);
	dec = dec*rad_per_deg;
	radectransform(ecinc, ra, dec, &eclon, &eclat);
	radectransform(invinc, (ra-invnode), dec, &invlon, &invlat);
	/* radecfancy(invnode, invinc, ra, dec, &invlat, &invlon); */
	ra = ra*deg_per_rad/deg_per_hrs;
	dec = dec*deg_per_rad;
	printf("%f\t%f\t%f\t%f\t%f\t%f\n", ra, dec, eclon, eclat, invlon, invlat);
      }
    }
	
    else printf("Unknown option\n");
  }
  else print_help();
  return(0);
}

void radectransform(float inc, float ra, float dec,
		    float *lon, float *lat) {
  /* given coords (ra&dec/lon&lat) in rad, returns lat/long (ra/dec) in deg */
  float temp;
  float x,y,z;
  float xp, yp, zp;
  
  x = cos(ra) * cos(dec);
  y = sin(ra) * cos(dec);
  z = sin(dec);
  
  xp = x; /* rotating around x = ra axis */
  yp = cos(inc)*y + sin(inc)*z;
  zp = -sin(inc)*y + cos(inc)*z;
  
  *lat = asin(zp)*deg_per_rad;
  /*  *lat = atan(zp/(sqrt(xp*xp+yp*yp)))*deg_per_rad; */
  /* lon needs to be put into right quadrant */
  /* xp is like cosine, yp is like sin of angle we're looking for */
  if (xp==0) {
    if (yp >0) 
      *lon = 90;
    else if (yp<0)
      *lon = 270;
    else *lon = 0 ; /*both x and y = 0 */
  }
  else
    *lon = atan(yp/xp)*deg_per_rad;
  if (xp < 0) /* out of range of atan function */
    *lon = *lon+180;
  if (*lon < 0 ) 
    *lon = *lon +360; /* just spin it around the dial once */
  return;
}

/* not used - just trying a slightly different way to do this - */
void radecfancy(float node, float inc, float ra, float dec,
		float *lon, float *lat) {
  float temp;
  float lonsin, loncos;
  
  temp = sin(node)*sin(inc)*cos(ra)*cos(dec)
    - cos(node)*sin(inc)*sin(ra)*cos(dec) + cos(inc)*sin(dec);
  *lat = asin(temp);
  
  lonsin = -cos(inc)*sin(node)*cos(ra)*cos(dec)
    +cos(node)*cos(inc)*sin(ra)*cos(dec) + sin(inc)*sin(dec);
  lonsin = lonsin / cos(*lat);
  
  loncos = cos(node)*cos(ra)*cos(dec)+sin(node)*sin(ra)*cos(dec);
  loncos = loncos / cos(*lat);
  
  *lat = *lat * deg_per_rad;
  
  if (lonsin >=0 ) {
    /* in first or second quadrant */
    *lon = acos(loncos);
    *lon = *lon * deg_per_rad;
  }
  else {
    if (loncos >=0 ) {
      /* in fourth quadrant */
      *lon = acos(loncos);
      *lon = *lon*deg_per_rad + 360;
    }
    else {
      /* third quadrant */
      *lon = asin(lonsin); /* this is negative of angle from left x */
      *lon = *lon*deg_per_rad;
      *lon = 180 - *lon;
    }
  }
  return;

}
