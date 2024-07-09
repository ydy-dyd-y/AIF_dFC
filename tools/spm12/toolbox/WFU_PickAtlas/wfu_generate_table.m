function retStat = wfu_generate_table(handles)

%--------------------------------------------------------------------------
% Function: wfu_generate_table(handles)
%
% Purpose:  This function is called when the user clicks the Generate Table
%           button in the wfu_pickatlas application. It does the following:
%           1. presents a dialog box to the user for selection of filenames.
%           2. If user makes a single selection, and it is a file with
%           extension '.flist', that file is assumed to contain a list of
%           Analyze .img files to read. These are then read and grouped into a
%           list.
%           3. In all other cases, the selected filenames are simply grouped into a list
%           4. The user is prompted to enter the name of an output ROI file
%           that will hold the table created by this function.
%           5. If the user did NOT select an .flist file, the file list is
%           saved as an .flist file with the same name as
%           that entered for the output ROI table.
%           6. The ROI table is created and saved.
%
% Calling Args:
%           handles = structure with handles and user data (see GUIDATA)
%
% Returns:  0 if processing completed successfully, 1 otherwise.
%
%==========================================================================
% C H A N G E   L O G
% BAS   06/04/04    Added check for existence of spm_get() and, if not
%                   found, calls wfu_pickfile() instead. Also added retStat
%                   (returned status indicator) so won't display 'Table
%                   Generated' message if cancel out of file selection
%                   dialog. Also replaced spm_input() with Matlab function inputdlg().
%--------------------------------------------------------------------------

global atlas_toolbox
retStat = 1; %in case get out early

PA = strcat(atlas_toolbox,	'/', handles.SelectedAtlasType.subdir, ...
					'/', handles.SelectedAtlasType.dispimage   );
%--------------------------------------------------------------------
% Display file selection dialog box and get user's choice(s):
%--------------------------------------------------------------------
    if (exist('spm_get'))
 	    P_list = spm_get( Inf, '*.img', 'Select list of Analyze image files' );
    else
        P_list = wfu_pickfile('*.img', 'Select list of Analyze image files');
    end
	if size( P_list, 1) == 0, 
		disp( 'No files selected. Table generation canceled.' );
		return; 
	end;

   
    %------------------------------------
    % If user gave a .flist file, use it:
    %------------------------------------
    if (size( P_list, 1)==1)
        [pathstr,name,ext] = fileparts(P_list(1,:));
        if (strcmp(strtok(ext),'.flist'))%note: strtok needed to remove any trailing blanks before test
            listOfFiles=P_list(1,:); %save the selected filename for the listOfFiles
            P_list = wfu_read_flist(listOfFiles);
        end
    end 
  
    %--------------------------------------------------------------
    % Display ROI table name dialog box to let user change it:
    %--------------------------------------------------------------
    ofid  = 1;
    default_name = strcat( 'ROI_Table_', datestr(now,30) );
    if (exist('spm_get'))
        set(handles.txtITD, 'String', 'Please select ROI table name');
        ofile_name = spm_input( 'Give a file name for the output: ', '1', 's', default_name);
        %Finter = spm_input('!GetWin');    %formerly spm_input_ui
        %delete(Finter);                    
        spm_figure('Clear','Interactive'); %clear input window
    else
        prompt = {'Give a file name for the output: '};
        dlg_title = 'Please select ROI table name';
        answer  = inputdlg(prompt,dlg_title,1,{default_name},'on');  
        if (isempty(answer))%if user canceled, return
            disp('User clicked cancel button. Table generation canceled.');
            return
        elseif (length(answer{1})==0)
            disp('User did not enter an output filename. Table generation canceled.');
            return        
        else
            ofile_name = answer{1};%else change from cell to array type for fopen call below.
        end
    end
    
    %---------------------------------
    % Open the output file for writing:
    %---------------------------------
    ofid = fopen( strcat( ofile_name, '.tbl'), 'w' );
    if (ofid==-1)
        errordlg(sprintf('Error opening file %s',strcat( ofile_name, '.tbl')));
        return
    end
    fprintf(ofid,'    Size\tAverage     \tStd.Dev.\tT        \t Region   \tROI name   \tL/R/B\tStudy     \tImage\n');
    disp(sprintf('Reading files, please wait...'));

    %------------------------------------------------
    % Create a .flist file from the ROI table name 
    % and write the list of filenames into it:
    %------------------------------------------------
    lfid  = 1;
    [pathstr,name] = fileparts(ofile_name);
    if (isempty(pathstr))
        pathstr=pwd;
    end
    listOfFiles = sprintf( '%s/%s.flist',pathstr,name);
    lfid = fopen(listOfFiles, 'w' );
    %------------------------------------------------
    % Write list (first record = number of filenames):
    %------------------------------------------------
    fprintf(lfid,'%d\n',size( P_list, 1));
    for ip = 1:size(P_list, 1)
        fprintf(lfid,'%s\n',P_list(ip,:));
    end
    fclose(lfid);
    
    %%%%%%%	Code borrowed from GenerateMask
	CMaskSide = [ 'R' 'L' 'B' ];
	List      = [];

	if (handles.isSimple)				% isSimple
	   AWL   = handles.WorkList;
	   AIndex      = 1:length( AWL );
	   for idx = AIndex
		A = AWL(idx).Atlas; 	R = AWL(idx).Region;	SubR = AWL(idx).Subregion;

		TheValue    = handles.Atlas(A).Region(R).SubregionValues(SubR);
		Offset      = handles.Atlas(A).Region(R).Offset;
		MaskSide    = handles.MaskSide;

		% Get the Point List and construct a Regn
		List    	= findindex( handles, TheValue, handles.Dilate, Offset, MaskSide);
    		Regn.names{1}	= handles.Atlas(A).Region(R).SubregionNames{SubR}; 
    		Regn.groups{1}	= handles.Atlas(A).Region(R).RegionName;

	    	print_ROI(ofid, List, Regn, CMaskSide(handles.MaskSide), PA, P_list, handles );
	   end
    	else % Not Simple = Advanced
           if(handles.display_flag==0) % NOT display_flag
             if(     ~isempty( handles.AdvancedWorkList ) &  ~isempty( handles.AdvancedIndex)   )
            	List = getadvancedlist(   handles.AdvancedWorkList( handles.AdvancedIndex(1)), handles);
            	Regn = getadvancedregion( handles.AdvancedWorkList( handles.AdvancedIndex(1)), handles);
	    	print_ROI(ofid, List, Regn, CMaskSide(handles.MaskSide), PA, P_list , handles);
             end
           else % display_flag
             if(handles.isAll)				% if isAll, then get all of finallist
                for i=1:length(handles.finallist)	
                    List = getadvancedlist(   handles.finallist(i).AdvancedWorkList, handles);
            	    Regn = getadvancedregion( handles.finallist(i).AdvancedWorkList, handles);
	    	    print_ROI(ofid, List, Regn, CMaskSide(handles.MaskSide), PA, P_list , handles);
             	end
             else 					%else if( ~isempty(CurrentFinal)) get that
                if(handles.CurrentFinal>0)
		  for i = 1:length(handles.CurrentFinal)	
                    List = getadvancedlist(   handles.finallist( handles.CurrentFinal(i) ).AdvancedWorkList, handles);
            	    Regn = getadvancedregion( handles.finallist( handles.CurrentFinal(i) ).AdvancedWorkList, handles);
	    	    print_ROI(ofid, List, Regn, CMaskSide(handles.MaskSide), PA, P_list , handles);
		  end
                end
             end % isAll
           end % is displayflag
	end % not isSimple
        
	fclose(ofid);
    retStat = 0;
	disp(sprintf(['Table written to ' pwd '/' ofile_name '.tbl']));
    return;
 

% --- Print ROI function
function print_ROI( ofid, reg_idx, Regn, Side, PA, P_list, handles )
% %%% 
	VA = spm_vol( PA );

	dim     = VA.dim;  plane = dim(1)*dim(2);
	reg_x   =      mod(reg_idx, dim(1));  % + 1 --> debugging found this off by one
	reg_y   =  fix(mod(reg_idx, plane ) / dim(1)) +1;
	reg_z   =  fix(    reg_idx/ plane ) +1;

	atlas_pix = [ reg_x, reg_y, reg_z, ones(length(reg_idx), 1) ]';
	atlas_mm  = VA.mat*atlas_pix; % VA.mat = pix2mm  

    nFiles = size(P_list, 1);
	for ip = 1:nFiles
        set(handles.txtITD, 'String', sprintf('processing file %d of %d...',ip,nFiles));
        drawnow; 
		PF = strtok(P_list( ip,:),' ');%strip trailing blanks
		VF = spm_vol( PF );
        %Get study ID from the pathname:
        [fpath fname fext] = fileparts( PF );
        [fstem fdir  fext] = fileparts( fpath );%back up one
		mm2pix   = inv( VF.mat);
		fmri_pix = mm2pix*atlas_mm;
        
		% hold = 1 --> trilinear interp; hold = 0 --> nearest neighbor
		% use 0 to debug when sampling original atlas, use 1 otherwise
		fmri_I   = spm_sample_vol( VF, fmri_pix(1,:), fmri_pix(2,:), fmri_pix(3,:), 1); 

		finite_idx = find( isfinite( fmri_I));
		if length( finite_idx) > 0,
			fmri_I = fmri_I( find( isfinite( fmri_I)));

			n_reg   = size(fmri_I,2); sum_reg = sum(fmri_I); ssq_reg = sum( fmri_I .* fmri_I);
			avg_reg = sum_reg/n_reg;  std_reg = sqrt(ssq_reg/n_reg - (avg_reg)^2);
			if std_reg > 0,         T_reg = avg_reg/std_reg;
			else,                   T_reg = sign(avg_reg)*Inf;      end;
            
			for ir = 1:length( Regn.groups )
				gRoups = union( Regn.groups(1), Regn.groups(ir) );
			end
			Region    = sprintf( '%s ', gRoups{:}      );
			Subregion = sprintf( '%s ', Regn.names{:}  );
		
			fprintf( ofid, ...
			      '%8g\t%8g\t%8g\t%8g\t%s\t%s\t%s\t%s\t%s\n',...
			       n_reg, avg_reg,std_reg,  T_reg, Region, Subregion, Side,   fdir, fname);
             
		end % fMRI_idx not empty
	end; % P_list 
%end print_ROI
    
