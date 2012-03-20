/* --------------------------------------------------- */
/*  Routine to read output netcdf file                 */
/*  analogue to the function IOIPSL restini+restget    */
/*  and to prepare initial condition for Tracer        */
/*  Model                                              */
/*  P.Lazzari & S.Salon - ECHO-OGA-INOGS October 2004  */
/* --------------------------------------------------- */
/* contains:                                           */
/* idx - function to inquire BC size                   */
/* idx_1d - function for 1d index vectors              */
/* bc_1d  - function for 1d BC    vectors              */
/* idx_3d - function for 3d index vectors              */
/* 2dt - function for 2d tracer variables              */
/* 3dt - function for 3d tracer variables              */
/* 2du - function for 2d bottom friction variables     */
/* 2dv - function for 2d bottom friction variables     */
/* 3du - function for 3d zonal      velocity variables */
/* 3dv - function for 3d meridional velocity variables */
/* 3dw - function for 3d vertical   velocity variables */
/* 3dFill - function to modify netcdf_files            */
/* rst - function to write restart netcdf_files        */
/* mesh_1D_f - function for meshmask 1d fields_Float   */
/* mesh_2D_f - function for meshmask 1d fields_Float   */
/* mesh_3D_f - function for meshmask 1d fields_Float   */
/* mesh_1D_d - function for meshmask 1d fields_Double  */
/* mesh_2D_d - function for meshmask 1d fields_Double  */
/* mesh_3D_d - function for meshmask 1d fields_Double  */
/* --------------------------------------------------- */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<netcdf.h>
#include<mpi.h>
/* Auxiliary function for Error Diagnostic */
void handle_error(char* nomevar ,int status);

/* ------------------------------------------------- */
/* function to inquire BC size */
int ioogsnc_idx(char *nomefile, char *nomedim, int *l ) 
{
      int i;
      int status, ncid, dimidp[1];
      int check=0;
      long int aux;
      double auxx;
      size_t lengthp[1];
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomedim,status);
       if (check) printf("IOOGS: bcfile= %s \n",nomefile);
       if (check) printf("IOOGS: dimension= %s \n",nomedim);
/*============================*/      

       status=nc_inq_dimid(ncid,nomedim,&dimidp[0]);
       status=nc_inq_dimlen(ncid,dimidp[0],lengthp);
       aux = (long int) lengthp[0];
       auxx = (double)  aux ;
       *l = (int) lengthp[0];
       handle_error(nomedim,status);
       if (check) printf("dim=%d \n",aux);
       if (check) printf("dim=%g \n",auxx);
       if (check) printf("dim=%d \n",*l);
       if (check) printf("dim=%d \n",lengthp[0]);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}
/* function for 1d integer variables */
int ioogsnc_idx_1d(char *nomefile, char *nomevar, int *ttpao, int *kkpao, int *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=1;
       count[3]=1;

       status=nc_get_vara_int(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %d \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */

/*============================*/      
return 0 ;
}


/* function for 1d integer variables */
int ioogsnc_idx_1d2(char *nomefile, char *nomevar, int *L, int *tem )
{
      int i;
      int status, ncid, dimidp[1];
      int  varidp,varidp1;
      size_t start[] = {0};
      size_t count[] = {0};
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/
       status=nc_open(nomefile,0,&ncid);
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=0;
       count[0]=*L;



       status=nc_get_vara_int(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %d \n",i, *(tem +i));
        }
       handle_error(nomevar,status);

/*=Close a NETcdf Dataset for Access: ncclose =======*/

       ncclose(ncid);
/* == free memory == */

/*============================*/
return 0 ;
}


/* function for 1d BC double variables */
int ioogsnc_bc_1d(char *nomefile, char *nomevar, int *ttpao, int *kkpao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=1;
       count[3]=1;

       status=nc_get_vara_double(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */

/*============================*/      
return 0 ;
}


/* function for 1d BC double variables */
int ioogsnc_bc_1d2(char *nomefile, char *nomevar, int *kkpao, double *tem )
{
      int i;
      int status, ncid; /*dimidp[4]; */
      int  varidp,varidp1;
      size_t start[] = {0};
      size_t count[] = {0};
/*      size_t lengthp[4];  */
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/
       status=nc_open(nomefile,0,&ncid);
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=0;
       count[0]=*kkpao;

       status=nc_get_vara_double(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);

/*=Close a NETcdf Dataset for Access: ncclose =======*/

       ncclose(ncid);
/* == free memory == */

/*============================*/
return 0 ;
}




/* function for 3d integer field */
int ioogsnc_idx_3d(char *nomefile, char *nomevar, int *ttpao, int *kkpao, int *jjpao, int *iipao, int *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;

       status=nc_get_vara_int(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %d \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 2d tracer variables */
int ioogsnc_2dt(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0}; 
      size_t count[] = {0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"deptht", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*jjpao;
       count[2]=*iipao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for reading 2d tracer variables */
int ioogsnc_2dt2(char *nomefile, char *nomevar, int *iipao, int *jjpao, float *tem )
{
      int i;
      int status, ncid, dimidp[2];
      int  varidp,varidp1;
      size_t start[] = {0,0};
      size_t count[] = {0,0};
      size_t lengthp[2];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/
       status=nc_open(nomefile,0,&ncid);
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/


       status=nc_inq_dimid(ncid,"y",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[0]);
/*============================*/

       status=nc_inq_dimid(ncid,"x",&dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[1]);
/*============================*/

    for(i=0;i<2;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=0;
       count[0]=*jjpao;
       count[1]=*iipao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);

/*=Close a NETcdf Dataset for Access: ncclose =======*/

       ncclose(ncid);

/*============================*/
return 0 ;
}



/* ------------------------------------------------- */
/* function for 3d tracer variables */
int ioogsnc_3dt(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, int *kkpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
     
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"deptht", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */

       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem);


    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 2d bottom friction variables */
int ioogsnc_2du(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0};
      size_t count[] = {0,0,0};
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;

/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"depthu", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      
    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */

       start[0]=(size_t)*ttpao-1;
       count[0]=1;
       count[1]=(size_t)*jjpao;
       count[2]=(size_t)*iipao;
       status=nc_get_vara_float(ncid, varidp, start, count, tem);

    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 2d bottom friction variables */
int ioogsnc_2dv(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0};
      size_t count[] = {0,0,0};
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;

     
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"depthv", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      
    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*jjpao;
       count[2]=*iipao;
       
       status=nc_get_vara_float(ncid, varidp, start, count, tem);

    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 3d zonal velocity variables */
int ioogsnc_3du(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, int *kkpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0};
      size_t count[] = {0,0,0,0};
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;

/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"deptht", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */

       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem);


    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 3d meridional velocity variables */
int ioogsnc_3dv(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, int *kkpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0};
      size_t count[] = {0,0,0,0};
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;

     
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"deptht", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;
       status=nc_get_vara_float(ncid, varidp, start, count, tem);

    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 3d vertical velocity variables */
int ioogsnc_3dw(char *nomefile, char *nomevar, int *ttpao, int *iipao, int *jjpao, int *kkpao, float *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0};
      size_t count[] = {0,0,0,0};
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
   
     
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time_counter",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"depthw", &dimidp[1]);
       if(status != NC_NOERR) {
       status=nc_inq_dimid(ncid,"deptht", &dimidp[1]);
       }	
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;


       status=nc_get_vara_float(ncid, varidp, start, count, tem);


    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      
       ncclose(ncid);
/*============================*/      
return 0 ;
}

/* ------------------------------------------------- */
/* function for 3d Tracer variables filling */
int ioogsnc_3dFill(char *nomefile, char *nomevar, int *iipao, int *jjpao, int *kkpao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t lengthp[4];
      static size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      nc_type xtypep;
      int check=0;
     
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,NC_WRITE,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*== Assign count Vector ==============*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
       count[i] = lengthp[i];
     }
    
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */      

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of tem data on the variable with name 'varidp' */

       status=nc_put_vara_double(ncid, varidp, start, count, tem);

    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);

/*============================*/      
return 0 ;
}

int ioogsnc(char *nomefile, char *nomevar, int *iipao, int *jjpao, int *kkpao, double *tem )
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t lengthp[4];
      nc_type xtypep;
      int check = 0;
      int rank;
      
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      
      
/*=Open a NETcdf Dataset for Access: nc_open =======*/
       if (check) printf("stage--->1");
       
       if (check) printf(":::initrc:::ioogsnc:::   rank No. %d opens file %s\n", rank, nomefile);
       status=nc_open(nomefile,0,&ncid);
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/
 
       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/
 
       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/
 
       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/
 
       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/
 
    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/
 
 
/*==Get a Variable ID from its name: nc_inq_varid*/
/* Tracer field*/
       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);
 
/*============================*/
 
       status=nc_get_var_double(ncid, varidp,tem);
 
    if(tem == NULL) printf("IOOGS:Error");
/*         for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
*/
       handle_error(nomevar,status);
/*=Close a NETcdf Dataset for Access: ncclose =======*/
 
       ncclose(ncid);
 
/*============================*/
return 0 ;
}
/* function for Tracer Restart Creation */

int ioogsnc_rst(char *nomefile, char *nmvarb, char *nmvarn, int *iipao, int *jjpao, int *kkpao, double *tem_lon, double *tem_lat,  double *tem_lev,   double *tem_tim, int * tstep, double *tem_info, double *tem_bck, double *tem_nxt )
{
      int i;
      int status, ncid, dimidp[7];
      int  varidp[8],varidp1;
      int dimlst1[1], dimlst2[2],dimlst3[3],dimlst4[4];
      size_t lengthp[7];
      size_t lat;
      size_t lon;
      size_t dep; 
      size_t tim;
      size_t tstp;
      size_t vardim[8];
      static size_t start1[] = {0};
      static size_t start2[] = {0,0};
      static size_t start3[] = {0,0,0};
      static size_t start4[] = {0,0,0,0};
      size_t count1[] = {0};
      size_t count2[] = {0,0};
      size_t count3[] = {0,0,0};
      size_t count4[] = {0,0,0,0};
      nc_type xtypep;
      
      int check = 0;      
/*=Open a NETcdf Dataset for Create: nc_create =======*/
       if (check) printf("stage--->1");
       status=nc_create(nomefile,0,&ncid);
       handle_error(nmvarb,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nmvarb);
/*=Dimension Definition=*/
       lengthp[0] = 1; 
       status=nc_def_dim(ncid, "time", lengthp[0], &dimidp[0]);
       handle_error(nmvarb,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/
 
       lengthp[1] = *kkpao; 
       status=nc_def_dim(ncid, "z", lengthp[1], &dimidp[1]);
       handle_error(nmvarb,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/
 
       lengthp[2] = *jjpao; 
       status=nc_def_dim(ncid, "y" ,lengthp[2], &dimidp[2]);
       handle_error(nmvarb,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/
 
       lengthp[3] = *iipao; 
       status=nc_def_dim(ncid, "x", lengthp[3], &dimidp[3]);
       handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/
 
       lengthp[4] = 1; 
       status=nc_def_dim(ncid, "x_a", lengthp[4], &dimidp[4]);
       handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[4]);
/*============================*/
 
       lengthp[5] = 1; 
       status=nc_def_dim(ncid, "y_a", lengthp[5], &dimidp[5]);
       handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[5]);
/*============================*/
 
       lengthp[6] = 3; 
       status=nc_def_dim(ncid, "z_a", lengthp[6], &dimidp[6]);
       handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[6]);
/*============================*/
 
    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nmvarb,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/
 
 
/*==Define  Variables : nc_def_dim */
/* Tracer field*/

       xtypep = NC_DOUBLE;
       vardim[0]=2;
       dimlst2[0]=dimidp[2];dimlst2[1]=dimidp[3];  
       status=nc_def_var(ncid,"nav_lon",xtypep,vardim[0],dimlst2,&varidp[0]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[0]);
 
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[1]=2;
       dimlst2[0]=dimidp[2];dimlst2[1]=dimidp[3];  
       status=nc_def_var(ncid,"nav_lat", xtypep ,vardim[1],dimlst2 ,&varidp[1]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[1]);
 
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[2]=1 ;
       dimlst1[0]=dimidp[1];
       status=nc_def_var(ncid,"nav_lev", xtypep, vardim[2], dimlst1, &varidp[2]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[2]);
 
/*============================*/
 
       xtypep = NC_DOUBLE;
       vardim[3]=1; 
       dimlst1[0]=dimidp[0];
       status=nc_def_var(ncid,"time", xtypep ,vardim[3], dimlst1 ,&varidp[3]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[3]);
 
/*============================*/
       
       xtypep = NC_INT;
       vardim[4]=1;
       dimlst1[0]=dimidp[0];
       status=nc_def_var(ncid, "time_step", xtypep , vardim[4] , dimlst1 , &varidp[4]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[4]);
 
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[5]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[6];
       dimlst4[2]=dimidp[5];
       dimlst4[3]=dimidp[4];
       status=nc_def_var(ncid, "info", xtypep ,vardim[5],dimlst4 ,&varidp[5]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[5]);
 
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[6]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvarb, xtypep ,vardim[6],dimlst4 , &varidp[6]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[6]);
 
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[7]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvarn, xtypep ,vardim[7],dimlst4 ,&varidp[7]);
       handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[7]);
 
/*============================*/
       nc_enddef(ncid); 
/*== Longitude ===============*/
/*
       count2[0]=lengthp[2];
       count2[1]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[0],start2,count2,tem_lon);
*/
       status=nc_put_var_double(ncid, varidp[0],tem_lon);

/*== Latitude ===============*/
       count2[0]=lengthp[2];
       count2[1]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[1],start2,count2,tem_lat);

/*== Level =================*/
       count1[0]=lengthp[1];
       status=nc_put_vara_double(ncid, varidp[2],start1,count1,tem_lev);

/*== Time =================*/
       count1[0]=lengthp[0];
       status=nc_put_vara_double(ncid, varidp[3],start1,count1,tem_tim);

/*== Time-Steps ===========*/
       count1[0]=lengthp[0];
       status=nc_put_vara_int(ncid, varidp[4],start1,count1,tstep);

/*== Info ================*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[4];
       count4[2]=lengthp[5];
       count4[3]=lengthp[6];
       status=nc_put_vara_double(ncid, varidp[5],start4,count4,tem_info);


/*== Trcbck ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[6],start4,count4,tem_bck);

/*== Trcnxt ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[7],start4,count4,tem_nxt);

    if(tem_nxt == NULL) printf("IOOGS:Error");
/*         for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nmvarb,status);
*/
/*=Close a NETcdf Dataset for Access: ncclose =======*/
 
       ncclose(ncid);
 
/*============================*/
return 0 ;
}


/* function for Tracer Restart Creation - New Name restart
 *  Elimino Time_step e Info
 *  Potrei aggiungere un attributo sul tempo, ma serve? */
int ioogsnc_rst2(char *nomefile, char *nmvarb, char *nmvarn, int *iipao, int *jjpao, int *kkpao, double *tem_lon, double *tem_lat,  double *tem_lev, double *tem_tim, double *tem_bck, double *tem_nxt )
{
      int i;
      int status, ncid, dimidp[7];
      int  varidp[8],varidp1;
      int dimlst1[1], dimlst2[2],dimlst3[3],dimlst4[4];
      static double miss_val_f[] = {1.e20};
      size_t lengthp[7];
      size_t lat;
      size_t lon;
      size_t dep;
      size_t tim;
      size_t tstp;
      size_t vardim[8];
      static size_t start1[] = {0};
      static size_t start2[] = {0,0};
      static size_t start3[] = {0,0,0};
      static size_t start4[] = {0,0,0,0};
      size_t count1[] = {0};
      size_t count2[] = {0,0};
      size_t count3[] = {0,0,0};
      size_t count4[] = {0,0,0,0};
      nc_type xtypep;
      int check = 0;
      char theTime[] = "19950115-00:00:00";
      strncpy(theTime, &nomefile[4], 17) ;
      printf("IOOGS THE_TIME: time is= %s and filename is %s \n",theTime, nomefile);

/*=Open a NETcdf Dataset for Create: nc_create =======*/
       if (check) printf("stage--->1");
       status=nc_create(nomefile,0,&ncid);
       handle_error(nmvarb,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nmvarb);
/*=Dimension Definition=*/
       lengthp[0] = 1;
       status=nc_def_dim(ncid, "time", lengthp[0], &dimidp[0]); handle_error(nmvarb,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/

       lengthp[1] = *kkpao;
       status=nc_def_dim(ncid, "z", lengthp[1], &dimidp[1]);handle_error(nmvarb,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/

       lengthp[2] = *jjpao;
       status=nc_def_dim(ncid, "y" ,lengthp[2], &dimidp[2]);handle_error(nmvarb,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/

       lengthp[3] = *iipao;
       status=nc_def_dim(ncid, "x", lengthp[3], &dimidp[3]); handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/

       lengthp[4] = 1;
       status=nc_def_dim(ncid, "x_a", lengthp[4], &dimidp[4]); handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[4]);
/*============================*/

       lengthp[5] = 1;
       status=nc_def_dim(ncid, "y_a", lengthp[5], &dimidp[5]); handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[5]);
/*============================*/

       lengthp[6] = 3;
       status=nc_def_dim(ncid, "z_a", lengthp[6], &dimidp[6]); handle_error(nmvarb,status);
       if (check) printf("XID=%d \n",dimidp[6]);
/*============================*/

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]); handle_error(nmvarb,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/


/*==Define  Variables : nc_def_dim */
/* Tracer field*/

       xtypep = NC_DOUBLE;
       vardim[0]=2;
       dimlst2[0]=dimidp[2];dimlst2[1]=dimidp[3];
       status=nc_def_var(ncid,"nav_lon",xtypep,vardim[0],dimlst2,&varidp[0]); handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[0]);

/*============================*/

       xtypep = NC_DOUBLE;
       vardim[1]=2;
       dimlst2[0]=dimidp[2];dimlst2[1]=dimidp[3];
       status=nc_def_var(ncid,"nav_lat", xtypep ,vardim[1],dimlst2 ,&varidp[1]); handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[1]);

/*============================*/

       xtypep = NC_DOUBLE;
       vardim[2]=1 ;
       dimlst1[0]=dimidp[1];
       status=nc_def_var(ncid,"nav_lev", xtypep, vardim[2], dimlst1, &varidp[2]); handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[2]);

/*============================*/

       xtypep = NC_DOUBLE;
       vardim[3]=1;
       dimlst1[0]=dimidp[0];
       status=nc_def_var(ncid,"time", xtypep ,vardim[3], dimlst1 ,&varidp[3]); handle_error(nmvarb,status);
       status=nc_put_att_text(ncid,varidp[3],"Units",33,"seconds since 1582-10-15 00:00:00");
       if (check) printf("variable Id=%d \n",varidp[3]);

/*============================*/
       xtypep = NC_DOUBLE;
       vardim[6]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvarb, xtypep ,vardim[6],dimlst4 , &varidp[6]); handle_error(nmvarb,status);
       status=nc_put_att_double(ncid,varidp[6] , "missing_value", xtypep, 1, miss_val_f); handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[6]);

/*============================*/

       xtypep = NC_DOUBLE;
       vardim[7]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvarn, xtypep ,vardim[7],dimlst4 ,&varidp[7]);  handle_error(nmvarb,status);
       status=nc_put_att_double(ncid,varidp[7] , "missing_value", xtypep, 1, miss_val_f); handle_error(nmvarb,status);
       if (check) printf("variable Id=%d \n",varidp[7]);

/*============================*/
       status=nc_put_att_text(ncid,NC_GLOBAL, "TimeString", 17 ,theTime); handle_error(nmvarb,status);

/*============================*/
       nc_enddef(ncid);
/*== Longitude ===============*/
/*
       count2[0]=lengthp[2];
       count2[1]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[0],start2,count2,tem_lon);
*/
       status=nc_put_var_double(ncid, varidp[0],tem_lon);

/*== Latitude ===============*/
       count2[0]=lengthp[2];
       count2[1]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[1],start2,count2,tem_lat);

/*== Level =================*/
       count1[0]=lengthp[1];
       status=nc_put_vara_double(ncid, varidp[2],start1,count1,tem_lev);

/*== Time =================*/
       count1[0]=lengthp[0];
       status=nc_put_vara_double(ncid, varidp[3],start1,count1,tem_tim);

/*== Time-Steps ===========*/
       count1[0]=lengthp[0];
/*       status=nc_put_vara_int(ncid, varidp[4],start1,count1,tstep); */


/*== Trcbck ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[6],start4,count4,tem_bck);

/*== Trcnxt ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];
       status=nc_put_vara_double(ncid, varidp[7],start4,count4,tem_nxt);

    if(tem_nxt == NULL) printf("IOOGS:Error");
/*         for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nmvarb,status);
*/
/*=Close a NETcdf Dataset for Access: ncclose =======*/

       ncclose(ncid);

/*============================*/
return 0 ;
}




int ioogsnc_ave_3d(char *nomefile, char *nmvar, int *iipao, int *jjpao, int *kkpao, float *tem_lon_f, float *tem_lat_f,  float *tem_lev_f,   double *tem_tim, float *tem_f)
{
      int i;
      int status, ncid, dimidp[7];
      static char convenctions[] = "COARDS";
      static float depth_range[] = {4.9991,4450.068};
      static float lat_range[] = {30.5,44.5};
      static float lon_range[] = {-9.25,36};
      static float miss_val_f[] = {1.e20};
      int  varidp[8],varidp1;
      int dimlst1[1], dimlst2[2],dimlst3[3],dimlst4[4];
      size_t lengthp[4];
      size_t lat;
      size_t lon;
      size_t dep; 
      size_t tim;
      size_t tstp;
      size_t vardim[8];
      static size_t start1[] = {0};
      static size_t start2[] = {0,0};
      static size_t start3[] = {0,0,0};
      static size_t start4[] = {0,0,0,0};
      size_t count1[] = {0};
      size_t count2[] = {0,0};
      size_t count3[] = {0,0,0};
      size_t count4[] = {0,0,0,0};
      nc_type xtypep;
/*      float tem_f[(*kkpao) * (*jjpao) *(*iipao)];
      float tem_lon_f[(*jjpao) *(*iipao)];
      float tem_lat_f[(*jjpao) *(*iipao)];
      float tem_lev_f[(*kkpao)];
      float *tem_f; 
      float *tem_lon_f;
      float *tem_lat_f;
      float *tem_lev_f;
*/
      int check = 1;      
/*============== Variables type conversion ===============*/
/*      tem_f = (float *)malloc (sizeof(float)* (*kkpao) * (*jjpao) *(*iipao));

        for(i=0;i<(*iipao)*(*jjpao)*(*kkpao);i++)
        {
         tem_f[i]= (float )*(tem +i) ;
        }
*/
/*============================*/
/*      tem_lon_f = (float *)malloc (sizeof(float)*(*iipao));

        for(i=0;i<(*iipao);i++)
        {
         tem_lon_f[i]= (float )*(tem_lon) + i* 0.125;
        }
*/
/*============================*/
/*      tem_lat_f = (float *)malloc (sizeof(float)*(*jjpao));
        for(i=0;i<(*jjpao);i++)
        {
         tem_lat_f[i]= (float )*(tem_lat) + i* 0.125;
        }
*/
/*============================*/
/*      tem_lev_f = (float *)malloc (sizeof(float)*(*kkpao));
        for(i=0;i<(*kkpao);i++)
        {
         tem_lev_f[i]= (float )*(tem_lev +i) ;
        }
*/
/*============================*/

/*=Open a NETcdf Dataset for Create: nc_create =======*/
       if (check) printf("stage--->1");
       status=nc_create(nomefile,0,&ncid);
       handle_error(nmvar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nmvar);
/*=Dimension Definition=*/
       lengthp[0] = 1; 
/*       status=nc_def_dim(ncid, "time", lengthp[0], &dimidp[0]);*/ 
       status=nc_def_dim(ncid, "time", NC_UNLIMITED, &dimidp[0]);
       handle_error(nmvar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/
 
       lengthp[1] = *kkpao; 
       status=nc_def_dim(ncid, "depth", lengthp[1], &dimidp[1]);
       handle_error(nmvar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/
 
       lengthp[2] = *jjpao; 
       status=nc_def_dim(ncid, "lat" ,lengthp[2], &dimidp[2]);
       handle_error(nmvar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/
 
       lengthp[3] = *iipao; 
       status=nc_def_dim(ncid, "lon", lengthp[3], &dimidp[3]);
       handle_error(nmvar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/
 
       lengthp[0] = 1; 
    for(i=1;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nmvar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/
 
 
/*==Define  Variables : nc_def_dim */
/* Tracer field*/
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[3]=1;
       dimlst1[0]=dimidp[0];
       status=nc_def_var(ncid,"time", xtypep ,vardim[3], dimlst1 ,&varidp[3]);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[3] , "units",7 ,"seconds");
       handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[3]);
/*============================*/

       xtypep = NC_FLOAT;
       vardim[2]=1 ;
       dimlst1[0]=dimidp[1];
       status=nc_def_var(ncid,"depth", xtypep, vardim[2], dimlst1, &varidp[2]);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[2] , "units",5 ,"meter");
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[2] , "positive",2 ,"up");
       handle_error(nmvar,status);
       status=nc_put_att_float(ncid, varidp[2] , "actual_range", xtypep, 2 ,depth_range);
       handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[2]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[1]=1;
       dimlst1[0]=dimidp[2];
       status=nc_def_var(ncid,"lat", xtypep ,vardim[1],dimlst1 ,&varidp[1]);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[1] , "units",13 ,"degrees_north");
       handle_error(nmvar,status);
       status=nc_put_att_float(ncid, varidp[1] , "actual_range", xtypep, 2, lat_range);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[1] , "long_name",8 ,"Latitude");
       handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[1]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[0]=1;
       dimlst1[0]=dimidp[3];
       status=nc_def_var(ncid,"lon",xtypep,vardim[0],dimlst1,&varidp[0]);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[0] , "units",12 ,"degrees_east");
       handle_error(nmvar,status);
       status=nc_put_att_float(ncid, varidp[0] , "actual_range", xtypep, 2, lon_range);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[0] , "long_name",9 ,"Longitude");
       handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[0]);
 
/*============================*/
       
       xtypep = NC_FLOAT;
       vardim[4]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvar, xtypep ,vardim[4],dimlst4 , &varidp[4]);
       handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[4] , "long_name",3 ,nmvar);
       handle_error(nmvar,status);
       status=nc_put_att_float(ncid, varidp[4] , "missing_value", xtypep, 1, miss_val_f);
       handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[6]);
 
/*== Global Attribute ========*/
       status=nc_put_att_text(ncid, NC_GLOBAL, "Convenctions",strlen(convenctions) ,convenctions);
       handle_error(nmvar,status);

       nc_enddef(ncid); 
/*== Longitude ===============*/

       count1[0]=lengthp[3];
       status=nc_put_vara_float(ncid, varidp[0],start1,count1,tem_lon_f);

/*== Latitude ===============*/
       count1[0]=lengthp[2];
       status=nc_put_vara_float(ncid, varidp[1],start1,count1,tem_lat_f);

/*== Level =================*/
       count1[0]=lengthp[1];
       status=nc_put_vara_float(ncid, varidp[2],start1,count1,tem_lev_f);

/*== Time =================*/
       count1[0]=lengthp[0];
       status=nc_put_vara_double(ncid, varidp[3],start1,count1,tem_tim);

/*== Trc ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];
/*        for(i=0;i<(*iipao)*(*jjpao)*(*kkpao);i++)
        {
         *(tem_f +i)= (float )*(tem +i) ;
        }
*/
       status=nc_put_vara_float(ncid, varidp[4],start4,count4,tem_f);

    if(tem_f == NULL) printf("IOOGS:Error");
/*         for(i=0;i<10;i++)
        {
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nmvarb,status);
*/
/*=Close a NETcdf Dataset for Access: ncclose =======*/
 
       ncclose(ncid);
 
/*============================*/
return 0 ;
}



int ioogsnc_ave_3d2(char *nomefile, char *nmvar, int *iipao, int *jjpao, int *kkpao, float *tem_lon_f, float *tem_lat_f,  float *tem_lev_f,   char *datestart, char* date__end, float *tem_f)
{
      int i;
      int status, ncid, dimidp[7];
      static char convenctions[] = "COARDS";
      static float depth_range[] = {4.9991,4450.068};
      static float lat_range[] = {30.5,44.5};
      static float lon_range[] = {-9.25,36};
      static float miss_val_f[] = {1.e20};
      int  varidp[8],varidp1;
      int dimlst1[1], dimlst2[2],dimlst3[3],dimlst4[4];
      size_t lengthp[4];
      size_t lat;
      size_t lon;
      size_t dep;
      size_t tim;
      size_t tstp;
      size_t vardim[8];
      static size_t start1[] = {0};
      static size_t start2[] = {0,0};
      static size_t start3[] = {0,0,0};
      static size_t start4[] = {0,0,0,0};
      size_t count1[] = {0};
      size_t count2[] = {0,0};
      size_t count3[] = {0,0,0};
      size_t count4[] = {0,0,0,0};
      nc_type xtypep;

      int check = 0;

/*=Open a NETcdf Dataset for Create: nc_create =======*/
       if (check) printf("stage--->1");
       status=nc_create(nomefile,0,&ncid);
       handle_error(nmvar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nmvar);
/*=Dimension Definition=*/
       lengthp[0] = 1;
       status=nc_def_dim(ncid, "time", NC_UNLIMITED, &dimidp[0]); handle_error(nmvar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/

       lengthp[1] = *kkpao;
       status=nc_def_dim(ncid, "depth", lengthp[1], &dimidp[1]); handle_error(nmvar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/

       lengthp[2] = *jjpao;
       status=nc_def_dim(ncid, "lat" ,lengthp[2], &dimidp[2]); handle_error(nmvar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/

       lengthp[3] = *iipao;
       status=nc_def_dim(ncid, "lon", lengthp[3], &dimidp[3]); handle_error(nmvar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/

       lengthp[0] = 1;
    for(i=1;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);handle_error(nmvar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/


/*==Define  Variables : nc_def_dim */
/* Tracer field*/
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[3]=1;
       dimlst1[0]=dimidp[0];
       status=nc_def_var     (ncid,"time", xtypep ,vardim[3], dimlst1 ,&varidp[3]);handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[3] , "units",7 ,"seconds");             handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[3]);
/*============================*/

       xtypep = NC_FLOAT;
       vardim[2]=1 ;
       dimlst1[0]=dimidp[1];
       status=nc_def_var(ncid,"depth", xtypep, vardim[2], dimlst1, &varidp[2]);handle_error(nmvar,status);
       status=nc_put_att_text (ncid, varidp[2] , "units",5 ,"meter"); handle_error(nmvar,status);
       status=nc_put_att_text (ncid, varidp[2] , "positive",2 ,"up");handle_error(nmvar,status);
       status=nc_put_att_float(ncid, varidp[2] , "actual_range", xtypep, 2 ,depth_range); handle_error(nmvar,status);

       if (check) printf("variable Id=%d \n",varidp[2]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[1]=1;
       dimlst1[0]=dimidp[2];
       status=nc_def_var(ncid,"lat", xtypep ,vardim[1],dimlst1 ,&varidp[1]); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[1], "units",13 ,"degrees_north"); handle_error(nmvar,status);
       status=nc_put_att_float(ncid,varidp[1], "actual_range", xtypep, 2, lat_range); handle_error(nmvar,status);
       status=nc_put_att_text( ncid,varidp[1], "long_name",8 ,"Latitude"); handle_error(nmvar,status);

       if (check) printf("variable Id=%d \n",varidp[1]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[0]=1;
       dimlst1[0]=dimidp[3];
       status=nc_def_var(ncid,"lon",xtypep,vardim[0],dimlst1,&varidp[0]); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[0] , "units",12 ,"degrees_east"); handle_error(nmvar,status);
       status=nc_put_att_float(ncid,varidp[0] , "actual_range", xtypep, 2, lon_range); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[0] , "long_name",9 ,"Longitude"); handle_error(nmvar,status);

       if (check) printf("variable Id=%d \n",varidp[0]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[4]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvar, xtypep ,vardim[4],dimlst4 , &varidp[4]); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[4] , "long_name",3 ,nmvar); handle_error(nmvar,status);
       status=nc_put_att_float(ncid,varidp[4] , "missing_value", xtypep, 1, miss_val_f); handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[6]);

/*== Global Attribute ========*/
       status=nc_put_att_text(ncid,NC_GLOBAL, "Convenctions",strlen(convenctions) ,convenctions); handle_error(nmvar,status);
       status=nc_put_att_text(ncid,NC_GLOBAL, "Time_Start", 17 ,datestart); handle_error(nmvar,status);
       status=nc_put_att_text(ncid,NC_GLOBAL, "Time___End", 17 ,date__end); handle_error(nmvar,status);

       nc_enddef(ncid);



/*== Longitude ===============*/

       count1[0]=lengthp[3];
       status=nc_put_vara_float(ncid, varidp[0],start1,count1,tem_lon_f);

/*== Latitude ===============*/
       count1[0]=lengthp[2];
       status=nc_put_vara_float(ncid, varidp[1],start1,count1,tem_lat_f);

/*== Level =================*/
       count1[0]=lengthp[1];
       status=nc_put_vara_float(ncid, varidp[2],start1,count1,tem_lev_f);

/*== Time ================= count1[0]=lengthp[0]; status=nc_put_vara_double(ncid, varidp[3],start1,count1,tem_tim);
*/


/*== Trc ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];

       status=nc_put_vara_float(ncid, varidp[4],start4,count4,tem_f);

    if(tem_f == NULL) printf("IOOGS:Error");

/*=Close a NETcdf Dataset for Access: ncclose =======*/

       ncclose(ncid);

/*============================*/
return 0 ;
}

/* writes ave backup files, temporary averages of variables */
int ioogsnc_ave_bkp(char *nomefile, char *nmvar, int *iipao, int *jjpao, int *kkpao, float *tem_lon_f, float *tem_lat_f,  float *tem_lev_f,   char *datestart, char* date__end, double *tem_f, int *avecounter)
{
      int i;
      int status, ncid, dimidp[7];
      static char convenctions[] = "COARDS";
      static float depth_range[] = {4.9991,4450.068};
      static float lat_range[] = {30.5,44.5};
      static float lon_range[] = {-9.25,36};
      static float miss_val_f[] = {1.e20};
      int  varidp[8],varidp1;
      int dimlst1[1], dimlst2[2],dimlst3[3],dimlst4[4];
      size_t lengthp[4];
      size_t lat;
      size_t lon;
      size_t dep;
      size_t tim;
      size_t tstp;
      size_t vardim[8];
      static size_t start1[] = {0};
      static size_t start2[] = {0,0};
      static size_t start3[] = {0,0,0};
      static size_t start4[] = {0,0,0,0};
      size_t count1[] = {0};
      size_t count2[] = {0,0};
      size_t count3[] = {0,0,0};
      size_t count4[] = {0,0,0,0};
      nc_type xtypep;

      int check = 0;

/*=Open a NETcdf Dataset for Create: nc_create =======*/
       if (check) printf("stage--->1");
       status=nc_create(nomefile,0,&ncid);
       handle_error(nmvar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nmvar);
/*=Dimension Definition=*/
       lengthp[0] = 1;
       status=nc_def_dim(ncid, "time", NC_UNLIMITED, &dimidp[0]); handle_error(nmvar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/

       lengthp[1] = *kkpao;
       status=nc_def_dim(ncid, "z", lengthp[1], &dimidp[1]); handle_error(nmvar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/

       lengthp[2] = *jjpao;
       status=nc_def_dim(ncid, "y" ,lengthp[2], &dimidp[2]); handle_error(nmvar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/

       lengthp[3] = *iipao;
       status=nc_def_dim(ncid, "x", lengthp[3], &dimidp[3]); handle_error(nmvar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/

       lengthp[0] = 1;
    for(i=1;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);handle_error(nmvar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/


/*==Define  Variables : nc_def_dim */
/* Tracer field*/
/*============================*/

       xtypep = NC_DOUBLE;
       vardim[3]=1;
       dimlst1[0]=dimidp[0];
       status=nc_def_var     (ncid,"time", xtypep ,vardim[3], dimlst1 ,&varidp[3]);handle_error(nmvar,status);
       status=nc_put_att_text(ncid, varidp[3] , "units",7 ,"seconds");             handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[3]);
/*============================*/

       xtypep = NC_FLOAT;
       vardim[2]=1 ;
       dimlst1[0]=dimidp[1];
       status=nc_def_var(ncid,"depth", xtypep, vardim[2], dimlst1, &varidp[2]);handle_error(nmvar,status);
       status=nc_put_att_text (ncid, varidp[2] , "units",5 ,"meter"); handle_error(nmvar,status);
       status=nc_put_att_text (ncid, varidp[2] , "positive",2 ,"up");handle_error(nmvar,status);
       status=nc_put_att_float(ncid, varidp[2] , "actual_range", xtypep, 2 ,depth_range); handle_error(nmvar,status);

       if (check) printf("variable Id=%d \n",varidp[2]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[1]=1;
       dimlst1[0]=dimidp[2];
       status=nc_def_var(ncid,"lat", xtypep ,vardim[1],dimlst1 ,&varidp[1]); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[1], "units",13 ,"degrees_north"); handle_error(nmvar,status);
       status=nc_put_att_float(ncid,varidp[1], "actual_range", xtypep, 2, lat_range); handle_error(nmvar,status);
       status=nc_put_att_text( ncid,varidp[1], "long_name",8 ,"Latitude"); handle_error(nmvar,status);

       if (check) printf("variable Id=%d \n",varidp[1]);

/*============================*/

       xtypep = NC_FLOAT;
       vardim[0]=1;
       dimlst1[0]=dimidp[3];
       status=nc_def_var(ncid,"lon",xtypep,vardim[0],dimlst1,&varidp[0]); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[0] , "units",12 ,"degrees_east"); handle_error(nmvar,status);
       status=nc_put_att_float(ncid,varidp[0] , "actual_range", xtypep, 2, lon_range); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[0] , "long_name",9 ,"Longitude"); handle_error(nmvar,status);

       if (check) printf("variable Id=%d \n",varidp[0]);

/*============================*/

       xtypep = NC_DOUBLE;
       vardim[4]=4;
       dimlst4[0]=dimidp[0];
       dimlst4[1]=dimidp[1];
       dimlst4[2]=dimidp[2];
       dimlst4[3]=dimidp[3];
       status=nc_def_var(ncid, nmvar, xtypep ,vardim[4],dimlst4 , &varidp[4]); handle_error(nmvar,status);
       status=nc_put_att_text (ncid,varidp[4] , "long_name",3 ,nmvar); handle_error(nmvar,status);
       status=nc_put_att_float(ncid,varidp[4] , "missing_value", xtypep, 1, miss_val_f); handle_error(nmvar,status);
       if (check) printf("variable Id=%d \n",varidp[6]);

/*== Global Attribute ========*/
       status=nc_put_att_text(ncid,NC_GLOBAL, "Convenctions",strlen(convenctions) ,convenctions); handle_error(nmvar,status);
       status=nc_put_att_text(ncid,NC_GLOBAL, "Time_Start", 17 ,datestart ); handle_error(nmvar,status);
       status=nc_put_att_text(ncid,NC_GLOBAL, "Time___End", 17 ,date__end ); handle_error(nmvar,status);
       status=nc_put_att_int (ncid,NC_GLOBAL, "Ave_counter", NC_INT, 1 ,avecounter); handle_error(nmvar,status);
       nc_enddef(ncid);



/*== Longitude ===============*/

       count1[0]=lengthp[3];
       status=nc_put_vara_float(ncid, varidp[0],start1,count1,tem_lon_f);

/*== Latitude ===============*/
       count1[0]=lengthp[2];
       status=nc_put_vara_float(ncid, varidp[1],start1,count1,tem_lat_f);

/*== Level =================*/
       count1[0]=lengthp[1];
       status=nc_put_vara_float(ncid, varidp[2],start1,count1,tem_lev_f);

/*== Time ================= count1[0]=lengthp[0]; status=nc_put_vara_double(ncid, varidp[3],start1,count1,tem_tim);
*/


/*== Trc ===============*/
       count4[0]=lengthp[0];
       count4[1]=lengthp[1];
       count4[2]=lengthp[2];
       count4[3]=lengthp[3];

       status=nc_put_vara_double(ncid, varidp[4],start4,count4,tem_f);

    if(tem_f == NULL) printf("IOOGS:Error");

/*=Close a NETcdf Dataset for Access: ncclose =======*/

       ncclose(ncid);

/*============================*/
return 0 ;
}




/* function for 1d meshmask float variables */
int ioogsnc_mesh_1d_f(char *nomefile, char *nomevar, int *ttpao, int *kkpao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0}; 
      size_t count[] = {0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
      float aux;
/*      float tem_f[*kkpao];*/
      float *tem_f;
      tem_f = (float *)malloc (sizeof(float)* (*kkpao) );

/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem_f);
       for(i=0;i<(*kkpao);i++)
         { aux = tem_f[i]; *(tem +i) = (double)  aux;}
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */
       free(tem_f);

/*============================*/      
return 0 ;
}
/* function for 2d meshmask float field */
int ioogsnc_mesh_2d_f(char *nomefile, char *nomevar, int *ttpao, int *jjpao, int *iipao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0}; 
      size_t count[] = {0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
      float aux;
/*      float tem_f[(*jjpao) * (*iipao)];*/
      float *tem_f;
      tem_f = (float *)malloc (sizeof(float) * (*jjpao) *(*iipao));

/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*jjpao;
       count[2]=*iipao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem_f);
       for(i=0;i<(*jjpao* *iipao);i++)
         { aux = tem_f[i]; *(tem +i) = (double)  aux;}
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */
       free(tem_f);
/*============================*/      
return 0 ;
}

/* function for 3d meshmask float field */
int ioogsnc_mesh_3d_f(char *nomefile, char *nomevar, int *ttpao, int *kkpao, int *jjpao, int *iipao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
      float aux;
/*      float tem_f[*kkpao * (*jjpao) * (*iipao)];*/
      float *tem_f;
      tem_f = (float *)malloc (sizeof(float)* (*kkpao) * (*jjpao) *(*iipao));

/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;

       status=nc_get_vara_float(ncid, varidp, start, count, tem_f);
       for(i=0;i<(*kkpao * *jjpao * *iipao);i++)
         { aux = tem_f[i]; *(tem +i) = (double)  aux;}
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */
       free(tem_f);
/*============================*/      
return 0 ;
}

/* function for 1d meshmask double variables */
int ioogsnc_mesh_1d_d(char *nomefile, char *nomevar, int *ttpao, int *kkpao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=1;
       count[3]=1;

       status=nc_get_vara_double(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */

/*============================*/      
return 0 ;
}
/* function for 2d meshmask double field */
int ioogsnc_mesh_2d_d(char *nomefile, char *nomevar, int *ttpao, int *jjpao, int *iipao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=1;
       count[2]=*jjpao;
       count[3]=*iipao;

       status=nc_get_vara_double(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */

/*============================*/      
return 0 ;
}

/* function for 3d meshmask double field */
int ioogsnc_mesh_3d_d(char *nomefile, char *nomevar, int *ttpao, int *kkpao, int *jjpao, int *iipao, double *tem ) 
{
      int i;
      int status, ncid, dimidp[4];
      int  varidp,varidp1;
      size_t start[] = {0,0,0,0}; 
      size_t count[] = {0,0,0,0}; 
      size_t lengthp[4];
      nc_type xtypep;
      int check=0;
/*=Open a NETcdf Dataset for Access: nc_open =======*/      
       status=nc_open(nomefile,0,&ncid);   
       handle_error(nomevar,status);
       if (check) printf("IOOGS: rstfile= %s \n",nomefile);
       if (check) printf("IOOGS: variabile= %s \n",nomevar);
/*============================*/      

       status=nc_inq_dimid(ncid,"time",&dimidp[0]);
       handle_error(nomevar,status);
       if (check) printf("TimeID=%d \n",dimidp[0]);
/*============================*/      

       status=nc_inq_dimid(ncid,"z", &dimidp[1]);
       handle_error(nomevar,status);
       if (check) printf("DepthID=%d \n",dimidp[1]);
/*============================*/      

       status=nc_inq_dimid(ncid,"y",&dimidp[2]);
       handle_error(nomevar,status);
       if (check) printf("YID=%d \n",dimidp[2]);
/*============================*/      

       status=nc_inq_dimid(ncid,"x",&dimidp[3]);
       handle_error(nomevar,status);
       if (check) printf("XID=%d \n",dimidp[3]);
/*============================*/      

    for(i=0;i<4;i++)
     {
       status=nc_inq_dimlen(ncid,dimidp[i],&lengthp[i]);
       handle_error(nomevar,status);
       if (check) printf("Dimension lenght[%d]=%d \n",dimidp[i],lengthp[i]);
     }
/*============================*/      
/* Get a Variable ID from its name: nc_inq_varid */

       status=nc_inq_varid(ncid,nomevar,&varidp);
       handle_error(nomevar,status);
       if (check) printf("variable Id=%d \n",varidp);

/* Transfer of read data with name 'varidp' on tem */
       start[0]=*ttpao-1;
       count[0]=1;
       count[1]=*kkpao;
       count[2]=*jjpao;
       count[3]=*iipao;

       status=nc_get_vara_double(ncid, varidp, start, count, tem);
    if(tem == NULL) printf("IOOGS:Error");
    if(check)
       for(i=0;i<10;i++)
        { 
         printf("IOOGS:tem[%d]= %g \n",i, *(tem +i));
        }
       handle_error(nomevar,status);
  
/*=Close a NETcdf Dataset for Access: ncclose =======*/      

       ncclose(ncid);
/* == free memory == */

/*============================*/      
return 0 ;
}


int get_att_int(char *nomefile,char*attname, int* res) {

    int status, ncid;
     status=nc_open(nomefile,0,&ncid); 	handle_error(nomefile,status);
     /*status=nc_get_att_text(ncid,NC_GLOBAL, "TimeString", 17 ,theTime); handle_error(nomefile,status);*/
     status=nc_get_att_int(ncid,NC_GLOBAL, attname, res); handle_error(nomefile,status);
     status=ncclose(ncid);handle_error(nomefile,status);
return 0 ;
}

int get_att_char(char *nomefile,char*attname, char* res) {

    int status, ncid;
     status=nc_open(nomefile,0,&ncid); 	handle_error(nomefile,status);
     status=nc_get_att_text(ncid,NC_GLOBAL, attname, res); handle_error(nomefile,status);
     status=ncclose(ncid);handle_error(nomefile,status);
return 0 ;
}





void handle_error(char *nomevar,int status)
{
  int error;
  if(status != NC_NOERR)
    {
      printf("IOOGS:%s \n",nomevar);      
      printf("%s \n",nc_strerror(status));
      exit(error);
    }
}
