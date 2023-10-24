/**
 * @file   hydralogon.C
 * @author Simon Lang
 * @date   Thu May 11 09:31:40 2006
 * 
 * @brief  hydralogin.C:
 *         A Root macro to make available all Hydra libraries in CINT.
 *         This macro can be listed in .rootrc
 * 
 * $Id: rootlogon.C,v 1.26 2008-05-21 13:54:48 halo Exp $
 *
 */


/** 
 * This function searches LD_LIBRARY_PATH to find each file of a 
 * list of dynamic libraries.
 * The default list of libraries is hard-coded here, and an additional 
 * list, or even a complete new list can be provided as arguments.
 * 
 * @param additional_libs All libraries in this space separated list are
 *        loaded after the default libraries.
 * @param load_from_second_location All libraries listed in this space 
 *        separated list are not loaded from the first location found in
 *        LD_LIBRARY_PATH, if they exist in another one, too.
 * @param personal_lib_list Overrides the default (hard-coded) list of
 *        libraries.
 * @param test_mode If this is true, then no libraries are loaded after
 *        the location and load order is shown.
 */
void rootlogon(TString additional_libs           = "FemtoCrap",
	       TString load_from_second_location = "",
	       TString personal_lib_list         = "",
	       Bool_t  test_mode                 = kFALSE)
{

   cout<<"WorkingDir "<<gSystem->WorkingDirectory()<<endl;

   TString common_libs;

   TString haddir = gSystem->GetFromPipe("echo $HADDIR");
   Bool_t useHADDIR =  ( haddir == "") ? kFALSE : kTRUE;

   TString myhaddir = gSystem->GetFromPipe("echo $MYHADDIR");
   Bool_t useMYHADDIR =  (myhaddir == "") ? kFALSE : kTRUE;

   if(!useHADDIR) cout<<"HADDIR is not set !"<<endl;


   // space seperated list and order of common libraries to be loaded
   common_libs  = "Hydra ";
    if((useHADDIR   && gSystem->AccessPathName(gSystem->ExpandPathName("${HADDIR}/lib/libRFIOtsm.so")) == 0 ) ||
       (useMYHADDIR && gSystem->AccessPathName(gSystem->ExpandPathName("${MYHADDIR}/lib/libRFIOtsm.so")) == 0)
     ) {
       common_libs += "RFIOtsm ";
   }
   common_libs += "Start ";
   common_libs += "Rich ";
   common_libs += "Mdc ";
   common_libs += "Tof ";
   common_libs += "Shower ";
   common_libs += "Emc ";
   common_libs += "Rpc ";
   common_libs += "Wall ";
   common_libs += "Sts ";
   common_libs += "FRpc ";
   common_libs += "Forward ";
   common_libs += "iTof ";
   common_libs += "PionTracker ";
   common_libs += "MdcTrackD ";
   common_libs += "MdcTrackG ";
   common_libs += "Kalman ";
   common_libs += "ShowerUtil ";
   common_libs += "Particle ";
   common_libs += "QA ";
   common_libs += "Dst ";
   common_libs += "Tools ";
   common_libs += "MdcUtil ";
   common_libs += "Simulation ";

   if((useHADDIR   && gSystem->AccessPathName(gSystem->ExpandPathName("${HADDIR}/lib/libOra.so")) == 0) ||
      (useMYHADDIR && gSystem->AccessPathName(gSystem->ExpandPathName("${MYHADDIR}/lib/libOra.so")) == 0)
     ) {
       common_libs += "Ora ";
       common_libs += "OraUtil ";
       common_libs += "OraSim ";
   }

   common_libs += "MdcGarfield ";
   common_libs += "Revt ";
   common_libs += "Online ";
   common_libs += "EventDisplay ";
   common_libs += "Alignment ";


   ////// here we go ...

   TString full_path;
   TString truncated_path;
   TString dir;
   TString file;
   TString filename;
   TString lib_name;
   TString path;
   TString path_dir;
   TString libs;
   TString lib;
   Bool_t  load_libs = kFALSE;
   Int_t   max_name  = 0;
   Bool_t  error     = kFALSE;

   //----------------------------------------------------------------
   // retrieve and clean DLL load path
   full_path = gSystem->GetDynamicPath();
   full_path.ReplaceAll( "::", ":" );
   full_path.Strip( TString::kBoth, ':' );
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // print path with one directory per line
   cout << "\nDLL Load Path Directories (in the Order of Precedence):" << endl;
   TObjArray* apath = full_path.Tokenize(":");
   if(apath){
       for(Int_t i = 0;i < apath->GetLast() + 1 ; i ++)
       {
	   path_dir = ((TObjString*) apath->At(i))->GetString();
           cout << "   " << path_dir.Data() << endl;
       }
       apath->Delete();
       delete apath;
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // handle personal list of libraries
   if (!personal_lib_list.IsNull())
   {
      common_libs = personal_lib_list;
   }
   common_libs.Prepend( " " );
   common_libs.Append( " " );
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // add a library from list of additional ones to common library list in
   // case it is not already a part of it - this check does not work if one
   // uses explicit paths for a library

   TString pluto_path = gSystem->Getenv( "PLUTODIR" );
   if (!pluto_path.IsNull()) {

       pluto_path = pluto_path + "/libPluto.so";

#if ROOT_VERSION_CODE  > ROOT_VERSION(6,0,0)
       pluto_path.ReplaceAll("//","/");
#else
       pluto_path.ReplaceAll("\/\/","\/");
#endif
       if(gSystem->AccessPathName(pluto_path.Data()) == 0) {

	   if(additional_libs.CompareTo("") == 0 ){ additional_libs = pluto_path; }
	   else                                   {
               additional_libs.Append( " " );
	       additional_libs += pluto_path;
	   }
       }
   }

   TObjArray* aaddlibs = additional_libs.Tokenize(" ");
   if(aaddlibs){
       for(Int_t i = 0;i < aaddlibs->GetLast() + 1 ; i ++)
       {
	   lib_name  = " ";
	   lib_name += ((TObjString*) aaddlibs->At(i))->GetString();
           lib_name += " ";

	   if (!common_libs.Contains( lib_name ))
	   {
	       common_libs += " ";
	       common_libs += lib_name;
	   }
       }
       aaddlibs->Delete();
       delete aaddlibs;
   }
   //----------------------------------------------------------------


   //----------------------------------------------------------------
   // find the longest filename
   TObjArray* acommonlibs = common_libs.Tokenize(" ");
   if(acommonlibs){
       for(Int_t i = 0;i < acommonlibs->GetLast() + 1 ; i ++)
       {
	   lib = ((TObjString*) acommonlibs->At(i))->GetString();
	   max_name = max_name < strlen( lib.Data() )  ?  strlen( lib.Data() )  :  max_name;
       }
       acommonlibs->Delete();
       delete acommonlibs;
   }
   //----------------------------------------------------------------

   load_from_second_location.Prepend( " " );
   load_from_second_location.Append( " " );
   cout << "\nList and Order of loaded DLLs and their Locations:" << endl;

   ProcessLibs:

   //----------------------------------------------------------------
   // get first library name

   acommonlibs = common_libs.Tokenize(" ");
   if(acommonlibs){
       for(Int_t i = 0;i < acommonlibs->GetLast() + 1 ; i ++)
       {
             lib = ((TObjString*) acommonlibs->At(i))->GetString();
	     filename = lib;
	     if (!filename.EndsWith( ".so" ))
	     {
		 filename  = "lib";
		 filename += lib;
		 filename += ".so";
	     }

	     //----------------------------------------------------------------
	     // determine full filename (of the first file location within path)
	     file = gSystem->Which( full_path.Data(), filename, kReadPermission );
	     if (file.IsNull())
	     {
		 if (!load_libs)
		 {
		     cerr << "  ERROR: " << filename << " not found!" << endl;
		     error = kTRUE;
		     continue;
		 }
		 return;
	     }
            //----------------------------------------------------------------

	     //----------------------------------------------------------------
	     // check, if the library should be loaded from a second path location
	     lib_name  = " ";
	     lib_name += lib;
	     lib_name += " ";
	     if (load_from_second_location.Contains( lib_name ))
	     {
		 // remove first location from path ...
		 dir             = ":";
		 dir            += gSystem->DirName( file.Data() );
		 dir            += ":";
		 truncated_path  = ":";
		 truncated_path += full_path;
		 truncated_path += ":";
		 truncated_path.ReplaceAll( dir, ":" );

		 // ... and search again, for a second location
		 filename =
		     gSystem->Which( truncated_path.Data(), filename, kReadPermission );
		 if (!filename.IsNull())
		 {
		     file = filename;
		 }
	     }
	     //----------------------------------------------------------------


	     //----------------------------------------------------------------
	     if (load_libs)
	     {
		 // load library ...
                 Int_t retval=0;
		 if( (retval = gSystem->Load( file.Data() )) != 0)
		 {
		     if(retval == -1 ){
			 cout << "ERROR: Failed to load " << lib  << " library! From : "<<file.Data()<<", Return code is "<<retval << endl;
			 return;
		     } else {
 			 cout << "INFO: Library " << lib  << " already in Memory!" << endl;

		     }
		 }
	     }
	     else
	     {
		 //----------------------------------------------------------------
		 // ... or just print the location
		 cout << "  ";
		 cout.setf( ios::left );
		 cout.width( max_name + 2 );
		 cout << lib_name << ": " << gSystem->DirName( file ) << endl;
		 //----------------------------------------------------------------
	     }
	     //----------------------------------------------------------------


       }
       acommonlibs->Delete();
       delete acommonlibs;
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // if everything is fine, process all libraries again, but now load them
   cout << endl;
   if (!load_libs  &&  !error  &&  !test_mode)
   {
      load_libs = kTRUE;
      goto ProcessLibs;
   }
   if ((!load_libs  &&  error)  ||  test_mode)
   {
      cout << "No libraries loaded!\n" << endl;
      return;
   }
   //----------------------------------------------------------------


   //----------------------------------------------------------------
   ////// ACLiC settings


   full_path = gSystem->Getenv( "MYHADDIR" );
   if (!full_path.IsNull())
   {
      full_path += "/include";
      if (!gSystem->AccessPathName( full_path.Data() ))
      {
	 full_path.Prepend( ".include " );
	 gROOT->ProcessLine( full_path.Data(), NULL );
         // gSystem->AddIncludePath( full_path.Data() ); somehow does not work!
      }
   }

   full_path = gSystem->Getenv( "HADDIR" );
   if (!full_path.IsNull())
   {
      full_path += "/include";
      if (!gSystem->AccessPathName( full_path.Data() ))
      {
	 full_path.Prepend( ".include " );
	 gROOT->ProcessLine( full_path.Data(), NULL );
         // gSystem->AddIncludePath( full_path.Data() ); somehow does not work!
      }
   }

   full_path = gSystem->Getenv( "PLUTODIR" );
   if (!full_path.IsNull())
   {
      full_path += "/include";
      if (!gSystem->AccessPathName( full_path.Data() ))
      {
	 full_path.Prepend( ".include " );
	 gROOT->ProcessLine( full_path.Data(), NULL );
         // gSystem->AddIncludePath( full_path.Data() ); somehow does not work!
      }
   }

   gSystem->SetFlagsOpt( "-O2" );

#if ROOT_VERSION_CODE  > ROOT_VERSION(6,0,0)
   #if ROOT_VERSION_CODE  > ROOT_VERSION(6,20,0)
      gSystem->SetMakeSharedLib(
        "cd $BuildDir; "
        "g++ -std=c++14 -c $Opt -pipe -Wall -fPIC -pthread $IncludePath $SourceFiles; "
        "g++ $ObjectFiles -shared -Wl,-soname,$LibName.so "
        "-O $LinkedLibs -o $SharedLib" );
   #else
      gSystem->SetMakeSharedLib(
         "cd $BuildDir; "
         "g++ -std=c++11 -c $Opt -pipe -Wall -fPIC -pthread $IncludePath $SourceFiles; "
         "g++ $ObjectFiles -shared -Wl,-soname,$LibName.so "
         "-O $LinkedLibs -o $SharedLib" );
   #endif
#else
   gSystem->SetMakeSharedLib(
      "cd $BuildDir; "
      "g++ -c $Opt -pipe -Wall -fPIC -pthread $IncludePath $SourceFiles; "
      "g++ $ObjectFiles -shared -Wl,-soname,$LibName.so "
      "-O $LinkedLibs -o $SharedLib" );
#endif



      //----------------------------------------------------------------


   //----------------------------------------------------------------
   ////// CINT/RINT settings


   // use spectral color palette for histograms
   gStyle->SetPalette( 1 );
   //----------------------------------------------------------------
}
