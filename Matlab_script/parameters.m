%% array definition

array = input(['Which array type was used for the recordings? \nChoose from ' ...
    'MidScala, SlimJ, CI612, CI622, CI632, Flex16, Flex20, Flex24, Flex28, ' ...
    'FlexSoft \nplace your response between single quotation marks: ']);

patient_number = input(['What is the number of your study participant, used for saving the result files? ' ...
    '\nplace your response between single quotation marks: ']);

%% insertion depth definition (in mm)

arr_EL_insdepth = {}; % insertion depth per electrode
arr_EL = {};

% Cochlear
nrEL = 22;

arr_EL.CI612 = 1:nrEL;
arr_EL_insdepth.CI612 = linspace(4.65, 18.9, nrEL); % active length = 14.25 mm, tip= 0.3 mm, distance from electrode tip to proximal white marker = 19.2 mm

arr_EL.CI622 = arr_EL.CI612;
arr_EL_insdepth.CI622 = linspace(5.6, 24.7, nrEL); % active length = 19.1 mm, tip= 0.3 mm, distance from electrode tip to proximal white marker = 25 mm

arr_EL.CI632 = arr_EL.CI612;
arr_EL_insdepth.CI632 = linspace(4, 18, nrEL); % active length = 14 mm, tip= 0.4 mm, distance from electrode tip to proximal white marker = 18.4 mm

% Advanced Bionics
nrEL = 16;

arr_EL.SlimJ = nrEL:-1:1;
arr_EL_insdepth.SlimJ = linspace(3, 22.5, nrEL); % distance from electrode tip to blue marker = 23 mm, contact spacing = 1.3 mm

arr_EL.MidScala = arr_EL.SlimJ;
arr_EL_insdepth.MidScala = linspace(3.375, 18, nrEL); % contact spacing = 0.975 mm

% MED-EL
nrEL = 12;

arr_EL.Flex16 = 1:nrEL;

arr_EL_insdepth.Flex16 = linspace(4, 16, nrEL); %contact spacing = 

arr_EL.Flex20 = arr_EL.Flex16;
arr_EL_insdepth.Flex20 = linspace(4.6, 20, nrEL); %active length = 15.4 mm --> contact spacing = 1.4 mm

arr_EL.Flex24 = arr_EL.Flex16;
arr_EL_insdepth.Flex24 = linspace(3.1, 24, nrEL); %active length = 20.9 mm --> contact spacing = 1.9 mm 

arr_EL.Flex28 = arr_EL.Flex16;
arr_EL_insdepth.Flex28 = linspace(4.9, 28, nrEL); %active length = 23.1 mm --> contact spacing = 2.1 mm

arr_EL.FlexSoft = arr_EL.Flex16;
arr_EL_insdepth.FlexSoft = linspace(5.1, 31.5, nrEL); %active length = 26.4 mm --> contact spacing = 2.4 mm

while ~isfield(arr_EL_insdepth, array)
    array = input(['\nThe chosen array type was not in the list of options \nChoose from ' ...
    'MidScala, SlimJ, CI612, CI622, CI632, Flex16, Flex20, Flex24, Flex28, FlexSoft' ...
    'FlexSoft \nplace your response between single quotation marks: ']);
end

insdepth = getfield(arr_EL_insdepth, array);
ELs = getfield(arr_EL, array);

%% insertion angle definition 

arr_EL_insangle = {};
nrEL = length(insdepth);

% Cochlear
arr_EL_insangle.CI612 = 390;
arr_EL_insangle.CI622 = 450;
arr_EL_insangle.CI632 = 390;

% Advanced Bionics
arr_EL_insangle.SlimJ = 425;
arr_EL_insangle.MidScala = 390;

% MED-EL
arr_EL_insangle.Flex16 = 270;
arr_EL_insangle.Flex20 = 360;
arr_EL_insangle.Flex24 = 450;
arr_EL_insangle.Flex28 = 580;
arr_EL_insangle.FlexSoft = 720;

ang = getfield(arr_EL_insangle, array);
insangle = linspace((insdepth(1)/insdepth(end))*ang, ang, nrEL);

%% cochlear duct length

axes_params = {};

axes_params.cochlear_length = input('\nWhat is the cochlear length in mm (press enter to use the default value)? ');

if isempty(axes_params.cochlear_length)
    axes_params.cochlear_length = 36.2;
    fprintf('\nThe cochlear length is set to the default value of %2.2f mm\n', axes_params.cochlear_length)
end

%% individualized axes

individualized_axes = input(['\nIf you want to use the default axes for this electrode type, ' ...
    'press 0. \nIf you want to define the axes using postoperative imaging, press 1. ']);

if individualized_axes

    axes_params.side = input(['\nWhat is the side of implantation?' ...
        '\nChoose between left and right' ...
        '\nplace your response between single quotation marks: ']);
    axes_params.image_name = input(['\nWhat is the image name, including extension?' ...
        '\nplace your response between single quotation marks: ']);
    axes_params.pixel_size = input('\nWhat is the pixel size in mm of your scan? ');
    % image dir - just place the image in the folder of the script
    % monitoring?

    while isempty(axes_params.image_name)
        axes_params.image_name = input(['\nImage name cannot be empty. What is the image name, including extension?' ...
            '\nplace your response between single quotation marks: ']);
    end
    
    while ~isfile(axes_params.image_name)
        axes_params.image_name = input(['\nImage does not exist in the directory. What is the image name, including extension?' ...
            '\nplace your response between single quotation marks: ']);    
    end
    
    while isempty(axes_params.pixel_size)
        axes_params.pixel_size = input('\nPixel size cannot be empty. What is the pixel size in mm of your scan? ');  
    end    

    while ~any([string('left'), 'right']==axes_params.side)
        axes_params.side = input(['\nSide of implantation should be left or right. ' ...
            '\nWhat is the side of implantation? ']);  
    end    

end
