#include "sys_info__q6_calc.h"

//This code modified from Rob's code
//
void sys_info__q6_calc::box_ti_in(box *i_box_ti) { box_ti = i_box_ti; }

void sys_info__q6_calc::create_bond_list() {

	int monomer_flag, number_particles;
	double bond_length_2;

	number_particles = box_ti->r_number_particles();

	bond_array.clear();
	bond_list.clear();

	bond_array.resize(number_particles);
	bond_list.resize(number_particles);

	for(int i = 0; i<number_particles; i++) {
		bond_array[i].resize(number_particles,0);
	}

	for(int part_1 = 0; part_1<number_particles; part_1++) {

		for(int part_2 = part_1 + 1; part_2<number_particles; part_2++) {

			if(box_ti->bond(part_1, part_2)) {

				bond_array[part_1][part_2] = 1;
				bond_array[part_2][part_1] = 1;

				bond_list[part_1].push_back(part_2);
				bond_list[part_2].push_back(part_1);

			}
		}
	}

	bond_list_exists = 1;
}

double sys_info__q6_calc::dxfn(int dim, int part_1, int part_2) {
	
	double delta;

	delta = box_ti->r_particle(part_1)->r_coordinate_box(dim) - box_ti->r_particle(part_2)->r_coordinate_box(dim);
	if(delta > box_ti->r_box_size(dim)*0.5) delta = delta - box_ti->r_box_size(dim);
	if(delta <= -box_ti->r_box_size(dim)*0.5) delta = delta + box_ti->r_box_size(dim);

	return delta;


}

void sys_info__q6_calc::q6_cor(double bin, double maxr, int q_r_flag, const char *tag) {

	int max = (int)(maxr/bin);
	int n, hist_all[max]; 
	double hist_q6[max];
	double d,dd,norm;
	int number_particles;

	number_particles = box_ti->r_number_particles();

	const static int ELL=6;
	double qimRe[number_particles][ELL+1], qimIm[number_particles][ELL+1];
	double cos_t, sin_t, phi;
	double rr, dx[3];
	int q6count,bondCount;
	int m;
	char diag=0;

	create_bond_list();

	// reset accumulators
	for(int i=0;i<number_particles;i++) for (m=0;m<=ELL;m++)
		qimRe[i][m]=qimIm[i][m]=0.0;

	// calculate and normalise the q for each m and each particle

	for(int i = 0; i<number_particles; i++) {
		for(int j = 0; j<bond_list[i].size(); j++) {

			// explicit EUCLIDEAN metric used here
			rr=0.0;
			for (int di=0;di<3;di++) {
				dx[di] = dxfn(di,i,bond_list[i][j]);
				rr += dx[di]*dx[di];
			}
			rr = sqrt(rr);
 
			cos_t = dx[2]/rr;
			sin_t = sqrt( 1 - (cos_t*cos_t) );
      
			if ( cos_t == 1.0 || cos_t == -1.0 ) phi=0.0;
			else
				phi = acos( dx[0] / sqrt( (dx[0]*dx[0]) + (dx[1]*dx[1] ) ) );
				// acos is in range[0..pi],  map to [-pi,pi]
				if ( dx[1] < 0.0 ) phi = -phi;
      
				if (diag) printf("c %.2e s %.2e p %.2e\n",cos_t,sin_t,phi);
      
				for (m=0;m<=ELL;m++) {
					rr = gsl_sf_legendre_sphPlm(ELL,m,cos_t);
					if (diag) printf("r%d %.2e\n",m,rr);
					qimRe[i][m] += rr * cos( m*phi );
					qimIm[i][m] += rr * sin( m*phi );
				}
		}

		// normalise, noting that sum from -m to +m is symmetric
		m=0;
		rr = (qimRe[i][m]*qimRe[i][m]) + (qimIm[i][m]*qimIm[i][m]) ;
		for (m=1;m<=ELL;m++) 
			rr += 2* ( (qimRe[i][m]*qimRe[i][m]) + (qimIm[i][m]*qimIm[i][m]) );
		rr = sqrt(rr);
		if (rr>0) 
			for (m=0;m<=ELL;m++) 
			{ qimRe[i][m]/=rr; qimIm[i][m]/=rr;  }

		if (diag) 
			for (m=0;m<=ELL;m++) 
				printf("qq %d %.2e %.2e\n",i,qimRe[i][m],qimIm[i][m]);

	}
  
	norm=0.0;  // this will be the total squared Q6
	for (m=0;m<=ELL;m++) {
		d=dd=0.0;
		for (int i=0;i<number_particles;i++) {
			d += qimRe[i][m];
			dd+= qimIm[i][m];
		}
		if (m==0) norm += (d*d) + (dd*dd);
		else      norm += 2.0 * ( (d*d) + (dd*dd) );
	}
	cout << box_ti->r_timestamp() << " " << norm << " " << number_particles << endl;

  for (int i=0;i<max;i++) { hist_all[i]=0; hist_q6[i]=0.0; }
  n=0;
  
  for (int i=0;i<number_particles;i++) for (int j=i+i;j<number_particles;j++) {
    d = box_ti->metric( i, j );
    if ( sqrt(d) <= maxr ) {
      d = sqrt(d)/bin;

      // dot prod of q6 of i and j
      m=0;
      rr = (qimRe[i][m]*qimRe[j][m]) + (qimIm[i][m]*qimIm[j][m]) ;
      for (m=1;m<=ELL;m++) {
        rr += 2* ( (qimRe[i][m]*qimRe[j][m]) + (qimIm[i][m]*qimIm[j][m]) );
      }

      hist_q6[ (int)d ] += rr;
      hist_all[ (int)d ] ++;

      n++;
    }
  }
 
  if(q_r_flag) {
	  for (int i=0;i<max;i++) {
	  	//a^3 - b^3 = (a-b)(a^2+ab+b^2) with a=i+1 and b=i
			  norm = (i+1)*(i+1) + i*(i+1) + i*i;
		  norm *= 0.5 * 4*M_PI/3.0*bin*bin*bin;
		  norm *= number_particles;
		  cout << (tag ? tag : "#q6 ") << (i+0.5)*bin << " " 
			  << hist_q6[i]/norm << " " 
			  << hist_all[i]/norm << " " 
			  << i << endl;
  	}	
  }

}
