%updated 4/3/2018 to include optional hand-drawn lesion mask. If left empty
%script will just use smri_dir & asl_base to perform PVc on ASL data
function ASLwarpPVc(real_lesion, smri_directory, asl_base)

disp(sprintf('Running ASL_RobustNormalization pipeline, date=%s...\n',datestr(now)))
disp(sprintf('real_lesion=%s...\n',real_lesion))
disp(sprintf('SMRIdir=%s...\n',smri_directory))
disp(sprintf('ASLdir=%s...\n',asl_base))

try
    for resource= dir(sprintf('%s/sequence*',asl_base))
        ASLoutput = sprintf('%s/%s',asl_base,resource.name);
        mkdir(fullfile(char(ASLoutput),'Normalized'));
        normfolder=fullfile(char(ASLoutput),'Normalized');
        %locate files in real_lesion
        if ~isempty(real_lesion)
            if isempty(spm_select('FPList',deblank(real_lesion),'^T1.nii$'))
                gunzip(spm_select('FPList',deblank(real_lesion),'^T1.nii.gz'));
            end
            realT1=spm_select('FPList',deblank(real_lesion),'^T1.nii$');
            
            if isempty(spm_select('FPList',deblank(real_lesion),'^lesionmask.nii$'))
                gunzip(spm_select('FPList',deblank(real_lesion),'^lesionmask.nii.gz'));
            end
            realles=spm_select('FPList',deblank(real_lesion),'^lesionmask.nii$');
        end
        
        %locate the files in ASLoutput
        CBF=spm_select('FPList',deblank(ASLoutput),'.*qCBF.*nii');CBF=CBF(1,:);
        V=spm_vol(CBF);nat_res=abs(diag(V.mat));nat_res=nat_res(1:3)';
        T1ASL=spm_select('FPList',deblank(ASLoutput),'^T1.nii');
        
        %locate desired files in RobustOUtput
        head=spm_select('FPList',smri_directory,'^head.nii$');%replace smri_directory with temporary folder in which the tar file contents are located
        
        %load matlabbatch file and perform
        %1) coregistration of T1ASL & T1lesmask to head.nii
        %2) Calculate GM/WM & nat2tpl field with VBM
        %3) warp CBF and lesionmask to tpl space
        %load('/projects/p20394/software/matlab/warp.mat');%replace with NUNDA location
        load('/projects/p20394/software/pipeline_external/ASL_RobustNormalization/coreg.mat');
      %  load('/home/yfc938/software/coreg.mat');
        matlabbatch{1,1}.spm.spatial.coreg.estwrite.ref{1}=head;
        matlabbatch{1,1}.spm.spatial.coreg.estwrite.source{1}=T1ASL;
        matlabbatch{1,1}.spm.spatial.coreg.estwrite.other=cellstr(CBF);
        if ~isempty(real_lesion)
            matlabbatch={matlabbatch{1,1},matlabbatch{1,1}};
            %matlabbatch{1,2}.spm.spatial.coreg.estimate.ref{1}=head;
            matlabbatch{1,2}.spm.spatial.coreg.estwrite.source{1}=realT1;
            matlabbatch{1,2}.spm.spatial.coreg.estwrite.other=cellstr(realles);
        end
        test2=load('/projects/p20394/software/pipeline_external/ASL_RobustNormalization/VBMest.mat');
        %test2=load('/home/yfc938/software/VBMest.mat');
        matlabbatch=cat(2,matlabbatch,test2.matlabbatch{1,1});
        matlabbatch{1,end}.spm.tools.vbm8.estwrite.data{1,1}=head;
        
        try
            spm_jobman('initcfg');
            spm_jobman('run_nogui',matlabbatch);
        end
        
        CBF=spm_select('FPList',deblank(ASLoutput),'^r.*qCBF.*nii');
        clear matlabbatch;
        matlabbatch=test2.matlabbatch(2);
        matlabbatch{1,1}.spm.tools.vbm8.tools.defs.field1=cellstr(spm_select('FPList',smri_directory,'^y_rhead.nii$'));
        matlabbatch{1,1}.spm.tools.vbm8.tools.defs.images=cellstr(CBF);
        if ~isempty(real_lesion)
            realT1=spm_select('FPList',deblank(real_lesion),'^rT1.nii$');
            matlabbatch{1,1}.spm.tools.vbm8.tools.defs.images=cat(1,cellstr(CBF),cellstr(realT1));
            matlabbatch={matlabbatch{1,1},test2.matlabbatch{1,3}};
            matlabbatch{1,2}.spm.tools.vbm8.tools.defs.field1=cellstr(spm_select('FPList',smri_directory,'^y_rhead.nii$'));
            realles=spm_select('FPList',real_lesion,'^rlesion.*nii');
            matlabbatch{1,2}.spm.tools.vbm8.tools.defs.images=cellstr(realles);
        end
        
        try
            spm_jobman('initcfg');
            spm_jobman('run_nogui',matlabbatch);
        end
        
        %Locate the warped template space tissue & brain masks
        GM=spm_read_vols(spm_vol(spm_select('FPList',smri_directory,'^wrp1.*.nii$')));
        WM=spm_read_vols(spm_vol(spm_select('FPList',smri_directory,'^wrp2.*.nii$')));
        CSF=spm_read_vols(spm_vol(spm_select('FPList',smri_directory,'^wrp3.*.nii$')));
        brain=GM+WM+CSF;brain(find(brain<0.5))=0;brain(find(brain))=1;
        brain(isnan(brain))=0;
        %located the template space ASL file and mask to >5
        ASLf=cellstr(spm_select('FPList',deblank(ASLoutput),'^wr.*qCBF.*nii$'));
        test=cellfun(@isempty,strfind(ASLf,'_PVc'));ASLf(find(test==0))=[];
        ASL=spm_read_vols(spm_vol(char(ASLf)));
        ASL(isnan(ASL))=0;
        ASL(find(ASL<5))=0;
        %Smooth tissue files to native resolution of ASL to match degree of
        %smoothness
        sGM=zeros(size(GM));spm_smooth(GM,sGM, nat_res);
        sWM=zeros(size(WM));spm_smooth(WM,sWM, nat_res);
        sWM99=sWM;sWM99(find(sWM99<0.99))=0;sWM99(find(sWM99))=1;
        %Limit calculations to GM>=0.3 to minimize elevated values due to division
        %by small numbers
        sGM(find(sGM<0.3))=0;
        V=spm_vol(spm_select('FPList',smri_directory,'^wrp1.*.nii$'));
        [pth,nm,ext]=fileparts(V.fname);
        V.fname=fullfile(normfolder,['wsGM',ext]);
        V.descrip='smoothed GM probability > 0.3';
        spm_create_vol(V);spm_write_vol(V,sGM);
        V.fname=fullfile(normfolder,['wsWM99',ext]);
        V.descrip='smoothed WM mask > 0.99';
        spm_create_vol(V);spm_write_vol(V,sWM99);
        
        %Perform partial volume correction by subtracting Pwm*WM_CBF and then
        %division by Pgm
        WMCBF=mean(ASL(find(sWM99)));
        ASL_PV=ASL-sWM*WMCBF;
        ASL_PV(find(sGM<0.3))=0;ASL_PV(find(sGM>=0.3))=ASL_PV(find(sGM>=.3))./sGM(find(sGM>=0.3));
        ASL_PV(isnan(ASL_PV))=0;
        if ~isempty(real_lesion)
            copyfile(spm_select('FPList',real_lesion,'^wrT1.nii$'),normfolder);
            V=spm_vol(spm_select('FPList',real_lesion,'wrlesionmask.nii$'));
            lesmask=spm_read_vols(V);
            lesmask(isnan(lesmask))=0;%lesmask(find(abs(lesmask)<0.01))=0;lesmask(find(lesmask))=1;
            lesmask(find(brain==0))=0;%lesmask=spm_erode(spm_dilate(lesmask));
            V.fname=fullfile(normfolder,'wrlesionmask.nii');V.descrip='lesion mask masked by brain>0.5 and SPMdilate&eroded';
            spm_create_vol(V);spm_write_vol(V,lesmask);
            brain(find(lesmask))=0;brain(find(CSF))=0;
            ASL_PV(find(lesmask))=0;
        end
        V=spm_vol(char(ASLf));
        [pth,nm,ext]=fileparts(V.fname);
        V.fname=fullfile(normfolder,[nm,'_PVc',ext]);V.descrip='Partial volume corrected (GM>=0.3) qCBF ';
        spm_create_vol(V);spm_write_vol(V,ASL_PV);
        
        
        %Define FOVmask, which combines (brain-lesion), GM>0.3 and ASL>5
        FOV=zeros(size(ASL));GMmask=sGM;GMmask(find(sGM))=1;
        FOV=GMmask.*brain;FOV(find(ASL<5))=0;
        V.dt=[4 0];V.descrip='tpl space FOVmask - 50% tissue & >30% smoothed GM';
        V.fname=fullfile(normfolder,['FOVmask',ext]);spm_create_vol(V);spm_write_vol(V,FOV);
        
        %copy files to ASLoutput
        copyfile(spm_select('FPList',smri_directory,'^wrp1.*.nii$'),normfolder);
        copyfile(spm_select('FPList',smri_directory,'^y_rhead.nii$'),normfolder);
        % ASL_PV=ASL;ASL_PV(find(sGM<0.3))=0;
        % ASL_PV(find(sGM>=0.3))=ASL_PV(find(sGM>=.3))./(sGM(find(sGM>=0.3))+0.4*sWM(find(sGM>=0.3)));
        % V.fname=fullfile(pth,[nm,'_PVcPET',ext]);V.descrip='Partial volume corrected (GM>=0.3) qCBF based on GM:WM=2.5';
    end
    
catch err
    disp(err.message);
    disp(err.identifier);
    for k=1:length(err.stack)
        fprintf('In %s at %d\n',err.stack(k).file,err.stack(k).line);
    end
    %exit;
end
return;
