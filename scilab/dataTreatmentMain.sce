clc
clear


function out = cycleGroups(toplevel,out)
 
    // consistency check
    if (~h5isFile(toplevel) & ~h5isGroup(toplevel)) then
        error("Wrong type of argument. HDF5 file or group expected.")
    end
    
    // cycles through groups with recurrecy
    for g = 1:size(toplevel.groups,1)
        currentObject = h5get(toplevel,toplevel.groups(g))
                
        disp("Analyzing... " + currentObject.path)
        
        if h5isGroup(currentObject) then
            out.name = currentObject.path
            out.groups(g) = cycleGroups(currentObject,out.groups(g))
        end
        
    end
    
endfunction




// main variables
path = "/home/leonardo/Projects/mbsim/2d-wagon-1.0/";
filename = "Wagon0.mbsim.h5";

cd(path);

// open H5 file
data = h5open(filename);

// create a tree of objects

// 1. create root node
root = uiCreateNode(data.name);
dummy = uiCreateNode("dummy");

// create the tree
h5tree = uiCreateTree(root,dummy);

// 2. create first level of nodes
out = struct()
out = cycleGroups(data.root,out)

h5tree = uiDeleteNode(h5tree,dummy);

uiDisplayTree(h5tree);

// close H5 file
h5close(data.root);
