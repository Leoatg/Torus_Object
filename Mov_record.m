
%Mov_record

classdef Mov_record < handle
    
    properties
        record_movie=0;
        record_mov_psi=0;
        record_mov_sxyz=0;
        re_psi_initialized=0;%external figure windows initialization
        re_sxyz_initialized=0;
        re_mov_begin;
        re_mov_end;
        re_mov_dt;
        
        GPU_calculation=0;
        Simulation_version='1st';
        
        currentFolder;
        dataFolder;
        mov_file_name_psi='Mov_psi_temp';
        mov_file_name_sxyz='Mov_sxyz_temp';
        
        vidObjpsi; F_psi;
        vidObjsxyz; F_sxyz;
        mov_profile='Motion JPEG AVI';
        Frame_Rate=4;
        zoom_factor=2;
        
        record_psi_initialized=0;%internal movie file initialization
        record_sxyz_initialized=0;
               
    end
    
    methods
        function Mov_psi(Movobj,fig_psi)
            if Movobj.record_movie==0
                fprintf('\nWarning: Movie record is not set.\n\n');
            else
                if Movobj.record_psi_initialized==0
                    if nargin<2
                        fprintf('\nWarning: No figue window. No Movie is recorded.\n\n');
                    else
                        set(fig_psi,'renderer','zbuffer');
                        figure(fig_psi);
                        if strcmp(Movobj.mov_file_name_psi,'Mov_psi_temp.avi')
                            fprintf('\nWarning: Movie name is not set.\n\n');
                        end
                        if isempty(Movobj.currentFolder) || isempty(Movobj.dataFolder)
                            fprintf('\nWarning: Movie saved path is not set.\n\n');
                        else
                            Movobj.mov_file_name_psi=strcat(Movobj.currentFolder,'\',Movobj.dataFolder,'\',Movobj.mov_file_name_psi);
                        end
                        Movobj.vidObjpsi=VideoWriter(Movobj.mov_file_name_psi,Movobj.mov_profile);
                        Movobj.vidObjpsi.FrameRate=Movobj.Frame_Rate;
                        open(Movobj.vidObjpsi);
                        Movobj.F_psi=[];
                        
                        Movobj.record_psi_initialized=1;
                    end
                    
                end
                if exist('fig_psi','var')
                    Movobj.F_psi = getframe(fig_psi);
                    writeVideo(Movobj.vidObjpsi,Movobj.F_psi);
                end
            end
        end %function Mov_psi end
        
        function Mov_sxyz(Movobj,fig_sxyz)
            if Movobj.record_movie==0
                fprintf('\nWarning: Movie record is not set.\n\n');
            else
                if Movobj.record_sxyz_initialized==0
                    if nargin<2
                        fprintf('\nWarning: No figue window. No Movie is recorded.\n\n');
                    else
                        set(fig_sxyz,'renderer','zbuffer');
                        figure(fig_sxyz);
                        if strcmp(Movobj.mov_file_name_sxyz,'Mov_sxyz_temp')
                            fprintf('\nWarning: Movie name is not set.\n\n');
                        end
                        if isempty(Movobj.currentFolder) || isempty(Movobj.dataFolder)
                            fprintf('\nWarning: Movie saved path is not set.\n\n');
                        else
                            Movobj.mov_file_name_sxyz=strcat(Movobj.currentFolder,'\',Movobj.dataFolder,'\',Movobj.mov_file_name_sxyz);
                        end
                        Movobj.vidObjsxyz=VideoWriter(Movobj.mov_file_name_sxyz,Movobj.mov_profile);
                        Movobj.vidObjsxyz.FrameRate=Movobj.Frame_Rate;
                        open(Movobj.vidObjsxyz);
                        Movobj.F_sxyz=[];
                        
                        Movobj.record_sxyz_initialized=1;
                    end
                    
                end
                if exist('fig_sxyz','var')
                    Movobj.F_sxyz = getframe(fig_sxyz);
                    writeVideo(Movobj.vidObjsxyz,Movobj.F_sxyz);
                end
            end
        end %function Mov_sxyz end
        
        function clear_all(Movobj)
            Movobj.re_psi_initialized=0;%for external figure windows
            Movobj.re_sxyz_initialized=0;
            Movobj.record_psi_initialized=0;%for internal movie file
            Movobj.record_sxyz_initialized=0;
            Movobj.vidObjpsi=[];
            Movobj.vidObjsxyz=[];
            Movobj.F_psi=[];Movobj.F_sxyz=[];
        end
    end
end