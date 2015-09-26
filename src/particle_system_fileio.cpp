/*
 * File: particle_system_fileio.cpp
 * Description: Implementation of particle system base class
 * for MC/MD sims.  File/IO functions for the class.
 * ---------------------------------------------------------
 *
 * John R. Dowdle
 * The University of Texas at Austin
 * June 2008
 * GPLv3.0+
 *
 */

#include "particle_system.h"



/******************* File I/O *******************************************/


void ParticleSystem::name_files(){
  coordfile_name = "";
  coordfile_name.append(out_dir);
  coordfile_name.append(proj_name);
  coordfile_name.append(".pdb");
  xdrcoordfile_name = "";
  xdrcoordfile_name.append(out_dir);
  xdrcoordfile_name.append(proj_name);
  xdrcoordfile_name.append(".xtc");
  corcoordfile_name = "";
  corcoordfile_name.append(out_dir);
  corcoordfile_name.append(proj_name);
  corcoordfile_name.append(".cor");
  tpifile_name = "";
  tpifile_name.append(out_dir);
  tpifile_name.append(proj_name);
  tpifile_name.append("_test.pdb");
  engfile_name = "";
  engfile_name.append(out_dir);
  engfile_name.append(proj_name);
  engfile_name.append(".eng");
  seqEngFileName = "";
  seqEngFileName.append(out_dir);
  seqEngFileName.append(proj_name);
  seqEngFileName.append(".seqEng");
  seqEngHistFileName = "";
  seqEngHistFileName.append(out_dir);
  seqEngHistFileName.append(proj_name);
  seqEngHistFileName.append(".seqEng.hist");
  seqCompHistFileName = "";
  seqCompHistFileName.append(out_dir);
  seqCompHistFileName.append(proj_name);
  seqCompHistFileName.append(".seqComp.hist");
  tpi_engfile_name = "";
  tpi_engfile_name.append(out_dir);
  tpi_engfile_name.append(proj_name);
  tpi_engfile_name.append(".weng");
  distfile_name = "";
  distfile_name.append(out_dir);
  distfile_name.append(proj_name);
  distfile_name.append(".cav.hist");
  seqEntropyHistFileName = "";
  seqEntropyHistFileName.append(out_dir);
  seqEntropyHistFileName.append(proj_name);
  seqEntropyHistFileName.append(".seqEnt.hist");
  bm_engfile_name = "";
  bm_engfile_name.append(out_dir);
  bm_engfile_name.append(proj_name);
  bm_engfile_name.append(".beng");
  growth_engfile_name = "";
  growth_engfile_name.append(out_dir);
  growth_engfile_name.append(proj_name);
  growth_engfile_name.append(".geng");
  logfile_name = "";
  logfile_name.append(log_dir);
  logfile_name.append(proj_name);
  logfile_name.append(".log");
  string rst_ext = ".rst";
  string rst_filename = proj_name + rst_ext;
  rst_path = rst_dir + rst_filename;
  cor_rst_name = rst_dir + proj_name + ".rst.cor";
  gr_cor_rst_name = rst_dir + proj_name + "_gr.rst.cor";

  if (rst_coords.length() > 0) {
    restart_coords = rst_dir + rst_coords;
  }
}


void ParticleSystem::assign_restart_coords(){
  vector<string> lines, line_data, tokens;
  string line;
  int anum, ptype;
  ifstream restart_file(restart_coords.c_str());
  map<int, int> ptypes_dummy;

  ptypes_dummy[0] = 0;  // LJ
  ptypes_dummy[1] = 0;  // HS
  ptypes_dummy[2] = 0;  // JG
  ptypes_dummy[3] = 0;  // SW

  if (restart_file.is_open()){
    while (!restart_file.eof()){
      getline(restart_file,line);
      Tokenize(line,tokens);
      if (tokens[0]=="FRAME") {
	//start_cycle = atoi(tokens[1].c_str());
	cubic_system = true;
	side_lengths.resize(ndims);
	for (int k=0; k<ndims; k++){
	  side_lengths[k] = atof(tokens[2+k].c_str());
	  if (k>0 && cubic_system){
	    if (side_lengths[k] != side_lengths[k-1]) cubic_system = false;
	  }
	}
	volume = get_system_volume();
	for (int i=0; i<natoms; i++){
	  getline(restart_file,line);
	  lines.push_back(line);
	}
	break;
      }
      tokens.clear();
    }
    if (lines.size()==0){
      fileio_error(restart_coords);
    } else {
      for (unsigned int i=0; i<lines.size(); i++){
	Tokenize(lines[i],line_data);
	anum = atoi(line_data[0].c_str());
	ptype = InverseParticleLabel(line_data[1]);
	ptypes_dummy[ptype]++;
	particles.push_back(Particle(ptype,hs_sigma,lj_sigma,lj_epsilon,
				     jg_lambda0, jg_lambda1, jg_lambda2,
				     jg_epsilon1, jg_epsilon2));
	for (unsigned int j=2; j<line_data.size(); j++){
	  gsl_matrix_set(coords,anum,j-2,atof(line_data[j].c_str()));
	}
	line_data.clear();
      }
    }
    restart_file.close();
  }
  else {
    fileio_error(restart_coords);
  }

  for (map<int,int>::iterator ii=ptypes_dummy.begin(); ii!=ptypes_dummy.end(); ++ii){
    ptypes[(*ii).first] = (*ii).second;
  }
}


void ParticleSystem::assign_bonds(){
  for (int i=0; i < polymer_size; i++){
    particles.push_back(Particle(1, hs_sigma));

    if (i == 0) {
      particles[i].bonds.push_back(1);
    } else if (i == polymer_size - 1) {
      particles[i].bonds.push_back(polymer_size - 2);
    } else {
      particles[i].bonds.push_back(i-1);
      particles[i].bonds.push_back(i+1);
    }
  }
}

void ParticleSystem::loadRstDOS(){
  ifstream dos_file(seqEntRstName.c_str());
  vector<string> tokens;
  string line;
  double engBin, value;

  if (dos_file.is_open()){
    while (!dos_file.eof()){
      getline(dos_file, line);
      Tokenize(line, tokens);
      if (tokens.size()==0){
	continue;
      } else if (tokens.size() >= 3) {
	if (tokens[0]=="#") {
	  if (tokens[1]=="gModFactor:"){
	    gModFactor = atof(tokens[2].c_str());
	  }
	} else {
	  engBin = (atof(tokens[0].c_str()) + atof(tokens[1].c_str()))/2.;
	  value = atof(tokens[2].c_str());
	  gsl_histogram_accumulate(seqEntropyHist, engBin, value);
	}
      }
      tokens.clear();
    }
  }  else {
    fileio_error(seqEntRstName);
  }
}


int ParticleSystem::ReadXTCFrame(XDRFILE *xdrtraj, string xdrtraj_path){
  int step, result;
  float time, precision;
  matrix box;
  rvec* xdrcoords;

  xdrcoords = (rvec*)calloc(natoms, sizeof(xdrcoords[0]));
  
  // open the file if it isn't already
  if (xdrtraj == NULL){
    // if a path is provided, try to open that path
    if (xdrtraj_path.length() > 0){
      xdrtraj = xdrfile_open(xdrtraj_path.c_str(),"r");
      // exit if we still can't open it
      if ( xdrtraj == NULL ){
	fprintf(logfile,"error:  unable to open the xdrfile\n");
	exit(FILEIO);
      }
    // if no path is provided, exit
    } else {
      fprintf(logfile,"error:  unable to open the xdrfile\n");
      exit(FILEIO);
    }
  }

  result = read_xtc(xdrtraj,natoms,&step,&time,box,xdrcoords,&precision);

  // if we haven't reached the end of the file
  if (result==0) {
    xdrbox2sidelengths(box);
    xdr2coords(xdrcoords);
  }

  volume = get_system_volume();
  frame_num = step;
  free(xdrcoords);

  return result;
}

// this is ugly.  don't blame me though.  blame
// VMD, catdcd, and FORTRAN
void ParticleSystem::ReadPDBFrame(ifstream &traj_file){
  vector<string> lines, line_data, tokens, atom_tokens;
  string line, atom_line;
  bool coords_loaded = false;
  int atoms_read=0;
  int pvs_frame_num=frame_num;
  //int cnt=0; // debug
  
  while(true){
    atoms_read = 0;
    if (coords_loaded || traj_file.eof()) break;
    getline(traj_file,line);
    Tokenize(line,tokens);         
    if (tokens.size()==0){
    } else if (tokens[0]=="CRYST1") {  // get the frame's box size
      for (int k=0; k<ndims; k++){
	side_lengths[k] = atof(tokens[1+k].c_str())/10.;
      }
      volume = get_system_volume();
    } else if (tokens[0]=="MODEL") { //<--we found a frame
      frame_num = atoi(tokens[1].c_str());
    } else if (tokens[0]=="ATOM") {  // load the current frame

      for (int k=0; k<ndims; k++){
	gsl_matrix_set(coords,atoms_read,k,atof(tokens[5+k].c_str())/10.);
      }

      atoms_read++;

      while(atoms_read < natoms){
	while (!traj_file.eof()){
	  getline(traj_file,atom_line);
	  Tokenize(atom_line,atom_tokens);
	  if (atom_tokens.size()>0){
	    if (atom_tokens[0]=="ATOM") break;
	  }
	}
	if (traj_file.eof()) break;
	for (int k=0; k<ndims; k++){
	  gsl_matrix_set(coords,atoms_read,k,atof(atom_tokens[5+k].c_str())/10.);
	}
	atom_tokens.clear();
	atoms_read++;
      }
      coords_loaded = true;
      if (pvs_frame_num==frame_num) frame_num++;
    }
    tokens.clear();

  }
}


void ParticleSystem::write_eng(){
  engfile = fopen(engfile_name.c_str(),"a");
  if (ncycles==0){
    string hdr_str0 = "# cycle  potential energy [kJ/mol]  ";
    string hdr_str1 = "volume [nm^3]  temperature [K]  ";
    string hdr_str2 = "pressure [kJ/mol/nm^-ndims]  virial [kJ/mol] \n";
    string hdr_str = hdr_str0 + hdr_str1 + hdr_str2;
    fprintf(engfile,"%s",hdr_str.c_str());
  }
  fprintf(engfile,"%d\t%.6e\t%0.6e\t%0.6e\t%.6e\t%.6e\n",	\
	  ncycles,pot_eng,volume,temp,pressure,virial);
  fclose(engfile);
}


void ParticleSystem::write_seq_eng(){
  string seqComp="";
  for (int i=0; i<seqLength; i++){
    if (particles[i].ptype==1)
      seqComp += "1";
    else
      seqComp += "0";
  }
  seqEngFile = fopen(seqEngFileName.c_str(),"a");
  if (ncycles==0){
    string hdr_str = "# cycle  seq. energy [kJ/mol]  \n";
    fprintf(seqEngFile,"%s", hdr_str.c_str());
  }
  fprintf(seqEngFile,"%d\t%.6e\t%s\n", ncycles, seqEng, seqComp.c_str());
  fclose(seqEngFile);
}


void ParticleSystem::talk(int fsnum){
  double part_acc_freq, part_zacc_freq, vol_acc_freq;

  fprintf(logfile,"# -------- System details ---------\n");
  fprintf(logfile,"# cycle: %d\n",(ncycles));
  fprintf(logfile,"# project name: %s\n",proj_name.c_str());
  fprintf(logfile,"# number of particles: %d\n",natoms);
  fprintf(logfile,"# particle types:  ");
  for (map<int,int>::iterator ii=ptypes.begin(); ii!=ptypes.end(); ++ii){
    fprintf(logfile,"%s",ParticleLabel((*ii).first).c_str());
    fprintf(logfile," (%d) ",(*ii).second);
  }
  fprintf(logfile,"\n");
  fprintf(logfile,"# number of dimensions:  %d\n",ndims);
  fprintf(logfile,"# potential eng:  %f\n",pot_eng);
  fprintf(logfile,"# temperature:  %0.2f\n",temp);
  if (ensemble=="npt") {
    fprintf(logfile,"# target pressure:  %0.2f\n",target_pressure); 
  }
  fprintf(logfile,"# current pressure:  %0.2f\n",pressure);
  fprintf(logfile,"# virial:  %0.2f\n",virial);
  fprintf(logfile,"# density:  %0.2f\n",density);
  fprintf(logfile,"# volume:  %0.2f\n",volume);
  for (int k=0; k<ndims; k++){
    fprintf(logfile,"# side length %d:  %0.2f\n",k,side_lengths[k]);
  }
  if (tail_cor){
    fprintf(logfile,"# tail correction (pres):  %0.2f\n",p_cor);
    fprintf(logfile,"# tail correction (eng):  %0.2f\n",eng_cor);
  }
  fprintf(logfile,"# cutoff:  %0.2f\n",cutoff);
  fprintf(logfile,"# number of cycles:  %d\n",(start_cycle+ncycles));
  //fprintf(logfile,"# number of accepted part. moves:  %d\n",npart_accepted);
  fprintf(logfile,"# delta_max:  %0.4f\n",delta_max);
  if (!cubic_system) fprintf(logfile,"# delta_zmax:  %0.4f\n",delta_zmax);
  if (ensemble=="npt"){
    if (cubic_system)
      fprintf(logfile,"# delta_lnV_max:  %0.4f\n",delta_lnV_max);
    else
      fprintf(logfile,"# delta_L_max:  %0.4f\n",delta_L_max);
  }
  if (ncycles > start_cycle){
    part_acc_freq = npart_accepted/(double)npart_attempted;
    fprintf(logfile,"# particle move acc. freq.:  %0.2f\n",part_acc_freq);
    if (!cubic_system){
      part_zacc_freq = npart_zaccepted/(double)npart_zattempted;
      fprintf(logfile,"# particle move zacc. freq.:  %0.2f\n",part_zacc_freq);
    }
    if (ensemble=="npt" && nvol_attempted > 0){
      vol_acc_freq = nvol_accepted/(double)nvol_attempted;
      fprintf(logfile,"# volume move acc. freq.:  %0.2f\n",vol_acc_freq);
    }
  }
  if (growth){
    //growth_sigma_avg = growth_sigma_sum / ng_att;
    fprintf(logfile,"# growth_sigma:  %0.2f\n",growth_sigma);
    //fprintf(logfile,"# growth_sigma_avg:  %0.2f\n",growth_sigma_avg);
    fprintf(logfile,"# p_growth:  %0.2f\n",p_growth);
  }
  fprintf(logfile,"# ---------------------------------\n");
}


void ParticleSystem::write_coords(){
  string pdb_line, cryst_line, aname;
  double rx, ry, rz;
  double slx, sly, slz;
  float precision=1000000;

  if (xdrfile) {
    xdrcoordfile = xdrfile_open(xdrcoordfile_name.c_str(),"a");
    matrix box;
    sidelengths2xdrbox(box);
    rvec *xdrcoords;
    xdrcoords = (rvec*)calloc(natoms, sizeof(xdrcoords[0]));
    coords2xdr(xdrcoords);
    write_xtc(xdrcoordfile, natoms, ncycles, ncycles, box,
	      xdrcoords, precision);
    xdrfile_close(xdrcoordfile);
    free(xdrcoords);
  } else if (corfile) {
    corcoordfile = fopen(corcoordfile_name.c_str(),"a");
    fprintf(corcoordfile,"FRAME %d\t",ncycles);
    for (int k=0; k<ndims; k++){
      fprintf(corcoordfile,"%0.6f\t",side_lengths[k]);
    }
    fprintf(corcoordfile,"\n");
    for (int i=0; i<natoms; i++) {
      fprintf(corcoordfile,"%d ",i); //particles[i].index);
      fprintf(corcoordfile," %s ",ParticleLabel(particles[i].ptype).c_str());
      for (int j=0; j<ndims-1; j++) {
	fprintf(corcoordfile,"%10.6f\t",gsl_matrix_get(coords,i,j));
      }
      fprintf(corcoordfile,"%10.6f\n",gsl_matrix_get(coords,i,ndims-1));
    }
    fprintf(corcoordfile,"ENDFRAME\n");
    fclose(corcoordfile);
  } else {
    coordfile = fopen(coordfile_name.c_str(),"a");
    fprintf(coordfile,"REMARK    GENERATED BY JLS\n");
    fprintf(coordfile,"TITLE     %s t=%10.2f\n",proj_name.c_str(),ncycles/1000.);
    fprintf(coordfile,"REMARK    THIS IS A SIMULATION BOX\n");
    cryst_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n";
    //cryst_line = "CRYST1   %6.3f   %6.3f   %6.3f  90.00  90.00  90.00 P 1           1\n";
    slx = side_lengths[0]*10.;
    sly = (ndims>1)?side_lengths[1]*10.:0.;
    slz = (ndims>2)?side_lengths[2]*10.:0.;
    fprintf(coordfile,cryst_line.c_str(),slx,sly,slz,90.0,90.0,90.0);
    fprintf(coordfile,"MODEL %8.1d\n",ncycles/cor_out_freq);
    pdb_line = "ATOM%7.1d  %s   %s  %4.1d    %8.3f%8.3f%8.3f  1.00  0.00      O\n";
    for (int i=0; i<natoms; i++){
      aname = ParticleLabel(particles[i].ptype).c_str();
      rx = gsl_matrix_get(coords,i,0)*10.;
      ry = (ndims>1)?gsl_matrix_get(coords,i,1)*10.:0.;
      rz = (ndims>2)?gsl_matrix_get(coords,i,2)*10.:0.;
      fprintf(coordfile,pdb_line.c_str(),i+1,aname.c_str(),aname.c_str(),i+1,rx,ry,rz);
    }
    fprintf(coordfile,"TER\n");
    fprintf(coordfile,"ENDMDL\n");
    fclose(coordfile);
  }
}


void ParticleSystem::write_tpi(){
  string pdb_line, cryst_line, aname;
  double rx, ry, rz;
  double slx, sly, slz;
  tpifile = fopen(tpifile_name.c_str(),"a");
  fprintf(tpifile,"REMARK    GENERATED BY JLS\n");
  fprintf(tpifile,"TITLE     %s t=%10.2f\n",proj_name.c_str(),ncycles/1000.);
  fprintf(tpifile,"REMARK    THIS IS A SIMULATION BOX\n");
  cryst_line = "CRYST1   %6.3f   %6.3f   %6.3f  90.00  90.00  90.00 P 1           1\n";
  slx = side_lengths[0]*10.;
  sly = (ndims>1)?side_lengths[1]*10.:0.;
  slz = (ndims>2)?side_lengths[2]*10.:0.;
  fprintf(tpifile,cryst_line.c_str(),slx,sly,slz);
  fprintf(tpifile,"MODEL %8.0d\n",ncycles/cor_out_freq);
  pdb_line = "ATOM%7.1d  %s   %s  %4.1d    %8.3f%8.3f%8.3f  1.00  0.00      O\n";
  for (int i=0; i<natoms; i++){
    aname = ParticleLabel(particles[i].ptype).c_str();
    rx = gsl_matrix_get(coords,i,0)*10.;
    ry = (ndims>1)?gsl_matrix_get(coords,i,1)*10.:0.;
    rz = (ndims>2)?gsl_matrix_get(coords,i,2)*10.:0.;
    fprintf(coordfile,pdb_line.c_str(),i+1,aname.c_str(),aname.c_str(),i+1,rx,ry,rz);
  }
  aname = ParticleLabel(particles[natoms].ptype).c_str();
  rx = gsl_matrix_get(coords,natoms,0)*10.;
  ry = (ndims>1)?gsl_matrix_get(coords,natoms,1)*10.:0.;
  rz = (ndims>2)?gsl_matrix_get(coords,natoms,2)*10.:0.;
  fprintf(tpifile,pdb_line.c_str(),natoms+1,aname.c_str(),aname.c_str(),natoms+1,rx,ry,rz);
  fprintf(tpifile,"TER\n");
  fprintf(tpifile,"ENDMDL\n");
  fclose(tpifile);
}


void ParticleSystem::write_tpi_eng(){
  engfile = fopen(tpi_engfile_name.c_str(),"a");
  if (ncycles==0){
    fprintf(engfile,"# frame  volume  nwsucc  nwtrials  <n>  <n^2>  pn\n");
  } else {
    if (calculation=="cavity" || calculation=="cavity_pmf"){
      fprintf(engfile,"%d\t%0.6e\t%d\t%d\t%0.6e\t%0.6e",frame_num,volume,
	      nwsucc,nwtrials,ncenters_mean,ncenters_sqrd_mean);
      for (int i=0; i<(int)pndist_cfg.size(); i++){
	fprintf(engfile,"\t%0.6e",pndist_cfg[i]);
      }
      fprintf(engfile,"\n");
    // } else if (calculation=="cavity_pmf") {
    //   double p0 = nwsucc/(double)nwtrials;
    //   double p0_cfg = nwsucc_cfg/(double)ninsertions;
    //   fprintf(engfile,"%d\t%0.6e\t%d\t%d\t%d\t%d\t%0.6e\t%0.6e\n",frame_num,volume,
    //          nwsucc_cfg,ninsertions,nwsucc,nwtrials,p0_cfg,p0);
    } else {
      fprintf(engfile,"%d\t%0.6e\t%d\t%d\t%0.6e\n",frame_num,volume,
	      nwsucc,nwtrials,tpi_eng_sum/double(ninsertions));
    }
  }
  fclose(engfile);
}


void ParticleSystem::write_cav_dist(){
  distfile = fopen(distfile_name.c_str(),"w");
  fprintf(distfile,"# %d\n",frame_num);
  gsl_histogram_fprintf(distfile, gsl_hist, "%g", "%g");
  fclose(distfile);
}

void ParticleSystem::write_seq_entropy_hist(){
  seqEntropyHistFile = fopen(seqEntropyHistFileName.c_str(),"w");
  fprintf(seqEntropyHistFile, "# %d\n", ncycles);
  fprintf(seqEntropyHistFile, "# gModFactor:  %g\n", gModFactor);
  gsl_histogram_fprintf(seqEntropyHistFile, seqEntropyHist, "%g", "%g");
  fclose(seqEntropyHistFile);
}

void ParticleSystem::write_seq_eng_hist(){
  seqEngHistFile = fopen(seqEngHistFileName.c_str(),"w");
  fprintf(seqEngHistFile, "# %d\n", ncycles);
  fprintf(seqEngHistFile, "# gModFactor:  %g\n", gModFactor);
  gsl_histogram_fprintf(seqEngHistFile, seqEngHist, "%g", "%g");
  fclose(seqEngHistFile);
}

void ParticleSystem::write_seq_comp_hist(){
  seqCompHistFile = fopen(seqCompHistFileName.c_str(),"w");
  fprintf(seqCompHistFile, "# %d\n", ncycles);
  for (int i=0; i < seqLength+1; i++){
    fprintf(seqCompHistFile, "# num HS:  %d\n", i);
    gsl_histogram_fprintf(seqCompHistFile, seqCompEngHists[i].hist, "%g", "%g");
  }
  fclose(seqCompHistFile);
}

void ParticleSystem::write_bm_eng(double bm_eng0, double bm_eng1){
  engfile = fopen(bm_engfile_name.c_str(),"a");
  if (ncycles==0){
    fprintf(engfile,"# frame  eng0  eng1\n");
  } else {
    fprintf(engfile,"%d\t%0.9e\t%0.9e\n",frame_num,bm_eng0,bm_eng1);
  }
  fclose(engfile);
}

void ParticleSystem::write_growth_eng(){
  engfile = fopen(growth_engfile_name.c_str(),"a");
  if (ncycles==0){
    fprintf(engfile,"# cycle  growth_sigma  ng_acc  ng_att  p_growth\n");
  } else {
    fprintf(engfile,"%d\t%0.6e\t%d\t%d\t%0.6e\n",ncycles,growth_sigma,
	    ng_acc,ng_att,p_growth);
  }
  fclose(engfile);
}


void ParticleSystem::write_rst(){
  rst_file = fopen(rst_path.c_str(),"w");
  //fprintf(rst_file,"proj_name = %s\n",proj_name);
  //fprintf(rst_file,"output_dir = %s\n",out_dir);
  //fprintf(rst_file,"ensemble = %s\n",ensemble);
  for (map<int,int>::iterator ii=ptypes.begin(); ii!=ptypes.end(); ++ii){
    fprintf(rst_file,"%s = %d\n",ParticleLabel((*ii).first).c_str(),\
	    (*ii).second);
  }
  fprintf(rst_file,"ndims = %d\n",ndims);
  fprintf(rst_file,"temp = %f\n",temp);
  if (ensemble=="npt"){
    fprintf(rst_file,"target_pressure = %f\n",target_pressure);
  }
  //fprintf(rst_file,"pressure = %f\n",pressure);
  //fprintf(rst_file,"virial = %f\n",virial);
  //fprintf(rst_file,"density = %f\n",density);
  //fprintf(rst_file,"volume = %f\n",volume);
  fprintf(rst_file,"cutoff = %f\n",cutoff);
  //if (tail_cor){
    //fprintf(rst_file,"p_cor =  %f\n",p_cor);
    //fprintf(rst_file,"eng_cor = %f\n",eng_cor);
  //}
  //fprintf(rst_file,"nmoves = %d\n",nmoves);
  fprintf(rst_file,"ncycles = %d\n",ncycles);
  //fprintf(rst_file,"npart_accepted = %d\n",npart_accepted);
  //fprintf(rst_file,"npart_attempted = %d\n",npart_attempted);
  fprintf(rst_file,"delta_max = %f\n",delta_max);
  if (!cubic_system) fprintf(rst_file,"delta_zmax = %f\n",delta_zmax);
  //fprintf(rst_file,"nvol_accepted = %d\n",nvol_accepted);
  //fprintf(rst_file,"nvol_attempted = %d\n",nvol_attempted);
  if (ensemble=="npt"){
    if (cubic_system)
      fprintf(rst_file,"delta_lnV_max = %f\n",delta_lnV_max);
    else
      fprintf(rst_file,"delta_L_max = %f\n",delta_L_max);
  }
  if (growth){
    //growth_sigma_avg = growth_sigma_sum / ng_att;
    fprintf(rst_file,"growth_sigma = %f\n",growth_sigma);
  }
  fclose(rst_file);
  write_rst_coords();
}


void ParticleSystem::write_rst_coords(){
  rst_coordfile = fopen(cor_rst_name.c_str(),"w");
  fprintf(rst_coordfile,"FRAME %d\t",ncycles);
  for (int k=0; k<ndims; k++){
    fprintf(rst_coordfile,"%0.6f\t",side_lengths[k]);
  }
  fprintf(rst_coordfile,"\n");
  for (int i=0; i<natoms; i++) {
    fprintf(rst_coordfile,"%d ",i); //particles[i].index);
    fprintf(rst_coordfile," %s ",ParticleLabel(particles[i].ptype).c_str());
    for (int j=0; j<ndims-1; j++) {
      fprintf(rst_coordfile,"%10.6f\t",gsl_matrix_get(coords,i,j));
    }
    fprintf(rst_coordfile,"%10.6f\n",gsl_matrix_get(coords,i,ndims-1));
  }
  fprintf(rst_coordfile,"ENDFRAME\n");
  fclose(rst_coordfile);
}


void ParticleSystem::write_growth_rst_coords(){
  rst_coordfile = fopen(gr_cor_rst_name.c_str(),"w");
  fprintf(rst_coordfile,"FRAME %d\t",ncycles);
  for (int k=0; k<ndims; k++){
    fprintf(rst_coordfile,"%f\t",side_lengths[k]);
  }
  fprintf(rst_coordfile,"\n");
  for (int i=0; i<natoms; i++) {
    fprintf(rst_coordfile,"%d ",i); //particles[i].index);
    fprintf(rst_coordfile," %s ",ParticleLabel(particles[i].ptype).c_str());
    for (int j=0; j<ndims-1; j++) {
      fprintf(rst_coordfile,"%10.4f\t",gsl_matrix_get(coords,i,j));
    }
    fprintf(rst_coordfile,"%10.4f\n",gsl_matrix_get(coords,i,ndims-1));
  }
  fprintf(rst_coordfile,"ENDFRAME\n");
  fclose(rst_coordfile);
}

