classdef Dataset
    
    properties
        fid = -1;
        filename = '';
        %datapath = '';
        grouppath = ''
        xmlpath = '';
        htypes = [];
    end
    
    methods
        
        function obj = Dataset(filename,groupname)
            
            % Set the hdf types
            obj.htypes = ismrmrd.util.hdf5_datatypes;
            
            % If the file exists, open it for read/write
            % otherwise, create it
            if exist(filename,'file')
                obj.fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            else
                fcpl = H5P.create('H5P_FILE_CREATE');
                obj.fid = H5F.create(filename,'H5F_ACC_TRUNC',fcpl,'H5P_DEFAULT');
                H5P.close(fcpl);
            end
            
            % Set the filename
            obj.filename = filename;
            
            % Set the group name
            %   default is dataset
            if nargin == 1
                groupname = 'dataset';
            end
            % Set the paths
            obj.grouppath = ['/' groupname];
            obj.xmlpath   = [obj.grouppath  '/xml'];
            %datapath  = [obj.grouppath '/data'];
            
            % Check if the group exists
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid,obj.grouppath,lapl_id) == 0)
                % group does not exist, create it
                group_id = H5G.create(obj.fid, obj.grouppath, 0);
                H5G.close(group_id);
            end
            H5P.close(lapl_id);
            
        end
        
        function obj = delete(obj)
            %Avoid corrupting the file due to data not written to disk
            obj.sync();
            % make sure we close the file when we clear or destroy the
            % class object
            obj.close();
        end
        
        function obj = close(obj)
            % close the file
            H5F.close(obj.fid);
        end
        
        function obj = sync(obj)
            % force writing all data to disk in HDF5 library cache,
            % hopefully avoids corruption
            H5F.flush(obj.fid,'H5F_SCOPE_GLOBAL');
        end
        
        function xmlstring = readxml(obj)
            % Check if the XML header exists
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid,obj.xmlpath,lapl_id) == 0)
                error('No XML header found.');
            end
            H5P.close(lapl_id);
            
            % Open
            xml_id = H5D.open(obj.fid, obj.xmlpath);
            
            % Get the type
            xml_dtype = H5D.get_type(xml_id);
            
            % Read the data
            hdr = H5D.read(xml_id, xml_dtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
            
            % Output depends on whether or not the stored string was variale length
            if (H5T.is_variable_str(xml_dtype))
                xmlstring = hdr{1};
            else
                xmlstring = hdr';
            end
            
            % Close the XML
            H5T.close(xml_dtype);
            H5D.close (xml_id);
        end
        
        function writexml(obj,xmlstring)
            % No validation is performed.  You're on your own.
            % make sure it's a char
            xmlstring = char(xmlstring);
            
            % TODO: add error checking on the write and return a status
            % TODO: if the matlab variable length string bug is resolved
            % then we should change this logic to just modify the length
            % and overwrite.
            
            % Check if the XML header exists
            %   if it does not exist, create it
            %   if it exists, modify the size appropriately
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid,obj.xmlpath,lapl_id) == 1)
                % Delete it
                H5L.delete(obj.fid, obj.xmlpath,'H5P_DEFAULT');
            end
            H5P.close(lapl_id);
            
            % Set variable length string type
            xml_dtype = H5T.copy('H5T_C_S1');
            % Matlab is having trouble writing variable length strings
            % that are longer that 512 characters.  Switched to fixed
            % length.
            H5T.set_size(xml_dtype,'H5T_VARIABLE');
            %H5T.set_size(xml_dtype, length(xmlstring));
            xml_space_id = H5S.create_simple (1, 1, []);
            xml_id = H5D.create(obj.fid, obj.xmlpath, xml_dtype, ....
                xml_space_id, 'H5P_DEFAULT');
            H5S.close(xml_space_id);
            
            % Write the data
            H5D.write(xml_id, xml_dtype, ...
                'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', {xmlstring});
            
            % Close the XML
            H5D.close(xml_id);
        end
        
        function nacq = getNumberOfAcquisitions(obj)
            datapath  = [obj.grouppath '/data'];
            % Check if the Data exists
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl_id) == 0)
                error([datapath ' does not exist in the HDF5 dataset.']);
            end
            dset = H5D.open(obj.fid, datapath);
            space = H5D.get_space(dset);
            H5S.get_simple_extent_dims(space);
            [~,dims,~] = H5S.get_simple_extent_dims(space);
            nacq = dims(1);
            H5S.close(space);
            H5D.close(dset);
            
        end
        
        
        
        function nacq = getNumberOfImages(obj, varname)
            if (~ischar(varname))
                error('varname must be a string');
            end
            % count the number of images by counting number of headers
            % based on ismrmrd_get_number_of_images()
            datapath  = [obj.grouppath '/' varname '/header'];
            % Check if the Data exists
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl_id) == 0)
                error([datapath ' does not exist in the HDF5 dataset.']);
            end
            dset = H5D.open(obj.fid, datapath);
            space = H5D.get_space(dset);
            H5S.get_simple_extent_dims(space);
            [~,dims,~] = H5S.get_simple_extent_dims(space);
            nacq = dims(1);
            H5S.close(space);
            H5D.close(dset);
            
        end
        
        function block = readAcquisition(obj, start, stop)
            if nargin == 1
                % Read all the acquisitions
                start = 1;
                stop = -1;
            elseif nargin == 2
                % Read a single acquisition
                stop = start;
            end
            datapath  = [obj.grouppath '/data'];
            % Check if the Data exists
            lapl=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl) == 0)
                error([datapath ' does not exist in the HDF5 dataset.']);
            end
            
            % Open the data
            dset = H5D.open(obj.fid, datapath);
            
            % Open the data space
            space = H5D.get_space(dset);
            
            % Get the size
            [~,dims,~] = H5S.get_simple_extent_dims(space);
            nacq = dims(1);
            
            % Create a mem_space for reading
            if (stop >= start)
                offset = [start-1];
                dims = [stop-start+1];
                mem_space = H5S.create_simple(1,dims,[]);
            else
                offset = [0];
                dims = [nacq];
                mem_space = H5S.create_simple(1,dims,[]);
            end
            
            % Read the desired acquisitions
            H5S.select_hyperslab(space,'H5S_SELECT_SET',offset,[1],[1],dims);
            d = H5D.read(dset, obj.htypes.T_Acquisition, ...
                mem_space, space, 'H5P_DEFAULT');
            
            % Pack'em
            block = ismrmrd.Acquisition(d.head, d.traj, d.data);
            
            % Clean up
            H5S.close(mem_space);
            H5S.close(space);
            H5D.close(dset);
        end
        
        function appendAcquisition(obj, acq)
            % Append an acquisition
            
            % TODO: Check the type of the input
            
            % The number of acquisitions that we are going to append
            N = acq.getNumber();
            
            % Create datapath string in acqusition functions rather than
            % class objects, to support ISMRMRD::Image as well
            datapath  = [obj.grouppath '/data'];
            % Check if the Data exists
            %   if it does not exist, create it
            %   if it does exist increase it's size
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl_id) == 0)
                % Data does not exist
                %   create with rank 1, unlimited, and set the chunk size
                dims    = [N];
                maxdims = [H5ML.get_constant_value('H5S_UNLIMITED')];
                file_space_id = H5S.create_simple(1, dims, maxdims);
                
                dcpl = H5P.create('H5P_DATASET_CREATE');
                chunk = [1];
                H5P.set_chunk (dcpl, chunk);
                data_id = H5D.create(obj.fid, datapath, ...
                    obj.htypes.T_Acquisition, ...
                    file_space_id, dcpl);
                H5P.close(dcpl);
                H5S.close(file_space_id);
                
            else
                % Open the data
                data_id = H5D.open(obj.fid, datapath);
                
                % Open the data space
                file_space_id = H5D.get_space(data_id);
                
                % Get the size, increment by N
                H5S.get_simple_extent_dims(file_space_id);
                [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
                dims = [dims(1)+N];
                H5D.set_extent (data_id, dims);
                H5S.close(file_space_id);
                
            end
            H5P.close(lapl_id);
            
            % Get the file space
            file_space_id = H5D.get_space(data_id);
            [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
            
            % Select the last N block
            offset = [dims(1)-N];
            H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',offset,[1],[1],[N]);
            
            % Mem space
            mem_space_id = H5S.create_simple(1,[N],[]);
            
            % Check and fix the acquisition header types
            acq.head.check();
            % TODO: Error checking on the sizes of the data and trajectories.
            
            % Pack the acquisition into the correct struct for writing
            d = struct();
            d.head = acq.head.toStruct();
            d.traj = acq.trajToFloat();
            d.data = acq.dataToFloat();
            
            % Write
            H5D.write(data_id, obj.htypes.T_Acquisition, ...
                mem_space_id, file_space_id, 'H5P_DEFAULT', d);
            
            % Clean up
            H5S.close(mem_space_id);
            H5S.close(file_space_id);
            H5D.close(data_id);
        end
        
        function attribString = readAttribString(obj,group)
            % Check if the attributes string exists
            grouppath = [obj.grouppath '/' group '/attributes'];
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid,grouppath,lapl_id) == 0)
                error('No attribute string found.');
            end
            H5P.close(lapl_id);
            
            % Open
            attrib_id = H5D.open(obj.fid, grouppath);
            
            % Get the type
            attrib_dtype = H5D.get_type(attrib_id);
            
            % Read the data
            hdr = H5D.read(attrib_id, attrib_dtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
            
            % Output depends on whether or not the stored string was variale length
            if (H5T.is_variable_str(attrib_dtype))
                attribString = hdr{1};
            else
                attribString = hdr';
            end
            
            % Close the attribute string
            H5T.close(attrib_dtype);
            H5D.close (attrib_id);
        end
        
         function appendAttribString(obj,group,attribString)
            % No validation is performed.  You're on your own.
            % make sure it's a char
            attribString = char(attribString);
            grouppath = [obj.grouppath '/' group '/attributes'];
            % TODO: add error checking on the write and return a status
            % TODO: if the matlab variable length string bug is resolved
            % then we should change this logic to just modify the length
            % and overwrite.
            
            % Check if the attribute string exists
            %   if it does not exist, create it
            %   if it exists, modify the size appropriately
            
            % TODO: Replace with code that can append attribute string
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid,grouppath,lapl_id) == 1)
                % Delete it
                H5L.delete(obj.fid, grouppath,'H5P_DEFAULT');
            end
            H5P.close(lapl_id);
            
            % Set variable length string type
            attrib_dtype = H5T.copy('H5T_C_S1');
            % Matlab is having trouble writing variable length strings
            % that are longer that 512 characters.  Switched to fixed
            % length.
            H5T.set_size(attrib_dtype,'H5T_VARIABLE');
            %H5T.set_size(attrib_dtype, length(attribString));
            attrib_space_id = H5S.create_simple (1, 1, []);
            attrib_id = H5D.create(obj.fid, grouppath, attrib_dtype, ....
                attrib_space_id, 'H5P_DEFAULT');
            H5S.close(attrib_space_id);
            
            % Write the data
            H5D.write(attrib_id, attrib_dtype, ...
                'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', {attribString});
            
            % Close the attribute
            H5D.close(attrib_id);
         end
        
        function imgHdr = readImageHeader(obj, varname,indx)

            datapath  = [obj.grouppath '/' varname '/header'];
            % Check if the Data exists
            lapl=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl) == 0)
                error([datapath ' does not exist in the HDF5 dataset.']);
            end
            
            % Open the data
            dset = H5D.open(obj.fid, datapath);
            
            % Open the data space
            space = H5D.get_space(dset);
            
            % Get the size
            [~,dims,~] = H5S.get_simple_extent_dims(space);
            nimg = dims(1);
            
            %nimg = obj.getNumberOfImages();
            % Check if index is within dimensions of the images set and
            % create a mem_space for reading
            %if (indx-1 < nimg )
            offset = [indx-1];
            dims = [nimg];
            mem_space = H5S.create_simple(1,1,dims);
            %else
            %    error('Requested Index exceeds number of images in dataset');
            %end
            % Read the desired acquisitions
            H5S.select_hyperslab(space,'H5S_SELECT_SET',offset,1,1,1);
            imgHdr = H5D.read(dset,obj.htypes.T_ImageHeader, ...
                mem_space, space, 'H5P_DEFAULT');
            
            
            % Clean up
            H5S.close(mem_space);
            H5S.close(space);
            H5D.close(dset);
        end
        
        function img = readImage(obj, varname, indx)
            %read ISMRMRD::Image from the dataset and variable specified in
            %varname
            if nargin < 2
                indx = 1;
            end
            %check varname input
            if (~ischar(varname))
                error('varname is not a string');
            end
            
            imgHdr = obj.readImageHeader(varname,indx);
            imgData =obj.readArray([varname '/data'],indx);
            imgAttribute = obj.readAttribString(varname);
            img = ismrmrd.Image(imgHdr,imgData,imgAttribute);
            
        end
        
        function appendImageHeader(obj,varname, hdr)
            %Append the ISMRMRD::ImageHeader
            % Create datapath string in acqusition functions rather than
            % class objects, to support ISMRMRD::Image as well
            datapath  = [obj.grouppath '/' varname '/header'];
            % Check if the Data exists
            %   if it does not exist, create it
            %   if it does exist increase it's size

            % We start with the header
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl_id) == 0)
                % Data does not exist
                %   create with rank 1, unlimited, and set the chunk size
                dims    = 1;
                maxdims = H5ML.get_constant_value('H5S_UNLIMITED');
                file_space_id = H5S.create_simple(1, dims, maxdims);
                
                dcpl = H5P.create('H5P_DATASET_CREATE');
                chunk = 1;
                H5P.set_chunk (dcpl, chunk);
                data_id = H5D.create(obj.fid, datapath, ...
                    obj.htypes.T_ImageHeader, ...
                    file_space_id, dcpl);
                H5P.close(dcpl);
                H5S.close(file_space_id);
                
            else
                % Open the data
                data_id = H5D.open(obj.fid, datapath);
                
                % Open the data space
                file_space_id = H5D.get_space(data_id);
                
                % Get the size, increment by N
                H5S.get_simple_extent_dims(file_space_id);
                [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
                dims = [dims(1)+1];
                H5D.set_extent (data_id, dims);
                H5S.close(file_space_id);
                
            end
            H5P.close(lapl_id);
            
            % Get the file space
            file_space_id = H5D.get_space(data_id);
            [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
            
            % Select the last N block
            offset = [dims(1)-1];
            H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',offset,[1],[1],1);
            
            % Mem space
            mem_space_id = H5S.create_simple(1,[1],[]);
            
            % Check and fix the acquisition header types
            hdr.check();
            d = hdr.toStruct();
            % TODO: Error checking on the sizes of the data and trajectories.
            
            % Pack the acquisition into the correct struct for writing
            
            
            % Write
            H5D.write(data_id,obj.htypes.T_ImageHeader, ...
                mem_space_id, file_space_id, 'H5P_DEFAULT', d);
            
            % Clean up
            H5S.close(mem_space_id);
            H5S.close(file_space_id);
            H5D.close(data_id);
        end
        
        function appendImage(obj, groupname,img)
            % Create an image
            
            % TODO: Check the type of the input
            
            % Create datapath string in acqusition functions rather than
            % class objects, to support ISMRMRD::Image as well
            grouppath = [obj.grouppath '/' groupname];
            
            % Create a new group to contain the new image
            plist = 'H5P_DEFAULT';
            g1 = H5G.open(obj.fid, obj.grouppath);
            if ~H5L.exists(g1,groupname,'H5P_DEFAULT');
                g1 = H5G.create(obj.fid, grouppath, plist,plist,plist);
            else
                error('Image already exists, TODO: add support for appending images');
            end
            
            datapath  = groupname;
            
            % Now we can append a new ImageHeader
            obj.appendImageHeader(groupname,img.head);
            
            % let's flush the HDF5 buffers to disk just in case
            obj.sync();
            
            % Convert Image Data to NDArray
            ndImage = ismrmrd.NDArray(img.data{1});
            
            % Using the high level interface for writing arrays
            
            %check varname input
            datapath = [grouppath '/data'];
            dims = size(ndImage.data);
            if length(dims) < 4
                dims(4) = 1; %Expand array to 4 columns
            elseif length(dims) > 4
                error('Image data has too many dimensions to be ISMRMRD::Image format');
            end
            
            [~,~,hltype] = ndImage.getNDArrayDataType();
            dims(5) = 1; %Expand to 5 dimensions
            dims(dims == 0) = 1; %Make sure that if we expanded to 5D, nothing is zero;
            
            [~,h5type,~] = ndImage.getNDArrayDataType();
            rank = length(dims);
            dims = fliplr(dims);
            start = zeros(size(dims));
            count = size(ndImage.data);
            count(rank) = 1;
            count = fliplr(count);
            sizeCount = count;
            
            sizeCount(1) = H5ML.get_constant_value('H5S_UNLIMITED');
            sizeCount(sizeCount == 0) = 1; 
            count(count == 0) = 1;
             %We want the object to be extendable along the  last (matlab) dimension, first hdf5 dimension
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl_id) == 0)
                % Data does not exist
                %   create with rank = numdimensions, unlimited, and set the chunk size
                
                file_space_id = H5S.create_simple(rank, dims, sizeCount);
                
                dcpl = H5P.create('H5P_DATASET_CREATE');
                chunk = count;
                H5P.set_chunk (dcpl, chunk);
                data_id = H5D.create(obj.fid, datapath, ...
                    h5type, ...
                    file_space_id, dcpl);
                H5P.close(dcpl);
                H5S.close(file_space_id);
                
            else
                % Open the data
                data_id = H5D.open(obj.fid, datapath);
                
                % Open the data space
                file_space_id = H5D.get_space(data_id);
                
                % Get the size, increment by N
                H5S.get_simple_extent_dims(file_space_id);
                [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
                dims(1) = dims(1)+1; %Appending a single NDArray
                H5D.set_extent (data_id, dims);
                H5S.close(file_space_id);
                start(1) = dims(1)-1;
            end
            H5P.close(lapl_id);
            
            % Get the file space
            file_space_id = H5D.get_space(data_id);
            [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
            %start(1) = dims(1);
            
            % Select the last N block
            %offset = start;
            H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',start,ones(rank,1),count,ones(rank,1));
            
            % Mem space
            mem_space_id = H5S.create_simple(rank,count,sizeCount);
           
            
            % Pack the acquisition into the correct struct for writing
            dataOut = ndImage.cplxToStruct();
            
            % Write
            H5D.write(data_id,h5type, ...
                mem_space_id, file_space_id, 'H5P_DEFAULT', dataOut);
            
            % Clean up
            H5S.close(mem_space_id);
            H5S.close(file_space_id);
            H5D.close(data_id);
            
            %append attribute
            obj.appendAttribString(groupname, img.attribute);
            
        end
        
        function data = readArray(obj, varname, indx)
            %read ISMRMRD::NDArray from the dataset and variable specified in
            %varname
            
            %check varname input
            if (~ischar(varname))
                error('varname is not a string');
            end
            
            %imgHdr = obj.readImageHeader(varname);
            datapath  = [obj.grouppath '/' varname];
            % Check if the Data exists
            lapl=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl) == 0)
                error([datapath ' does not exist in the HDF5 dataset.']);
            end
            
            % Open the data
            dset = H5D.open(obj.fid, datapath);
            
            % Get data type
            hdf5type = H5D.get_type(dset);
            
            % Open the data space
            space = H5D.get_space(dset);
            
            % Get the size
            [rankDims,dims,~] = H5S.get_simple_extent_dims(space);
            nArray = dims(1);
            count = dims;
            count(1) = 1;
            %nimg = obj.getNumberOfImages();
            % Check if index is within dimensions of the array set and
            % create a mem_space for reading
            if (indx-1 < nArray )
                %offset = [indx - 1];
                start = zeros(size(dims));
                start(1) = indx - 1;
                mem_space = H5S.create_simple(rankDims,count,[]);
            else
                error('Requested Index exceeds number of images in dataset');
            end
            % Read the desired acquisitions
            H5S.select_hyperslab(space,'H5S_SELECT_SET',start,[],[],count);
            dataTemp = H5D.read(dset, hdf5type, ...
                mem_space, space, 'H5P_DEFAULT');
            data = ismrmrd.NDArray(dataTemp);
            
            
            % Clean up
            H5S.close(mem_space);
            H5S.close(space);
            H5D.close(dset);
            
            
        end
        function appendArray(obj, varname, data)
            % write ISMRMRD::NDArray to the dataset and variable specified in
            % varname
            
            % Using the high level interface for writing arrays
            
            %check varname input
            if (~ischar(varname))
                error('varname is not a string');
            end
            
            if (~isa(data, 'ismrmrd.NDArray'))
                error('data is not an object of class ismrmrd.NDArray');
            end
            
            datapath = [obj.grouppath '/' varname];
            dims = size(data.data);
            [~,h5type,~] = data.getNDArrayDataType();
            dims(end+1) = 1; %Need to make sure the rank of the hdf5 dataspace is one greater than the
            %dimension of the array for appending
            dims = fliplr(dims);
            start = zeros(size(dims));
            count = size(data.data);
            count(end+1) = 1;
            count = fliplr(count);
            sizeCount = count;
            rank = length(dims);
            sizeCount(1) = H5ML.get_constant_value('H5S_UNLIMITED');
             %We want the object to be extendable along the  last (matlab) dimension, first hdf5 dimension
            lapl_id=H5P.create('H5P_LINK_ACCESS');
            if (H5L.exists(obj.fid, datapath, lapl_id) == 0)
                % Data does not exist
                %   create with rank = numdimensions, unlimited, and set the chunk size
                
                file_space_id = H5S.create_simple(rank, dims, sizeCount);
                
                dcpl = H5P.create('H5P_DATASET_CREATE');
                chunk = count;
                H5P.set_chunk (dcpl, chunk);
                data_id = H5D.create(obj.fid, datapath, ...
                    h5type, ...
                    file_space_id, dcpl);
                H5P.close(dcpl);
                H5S.close(file_space_id);
                
            else
                % Open the data
                data_id = H5D.open(obj.fid, datapath);
                
                % Open the data space
                file_space_id = H5D.get_space(data_id);
                
                % Get the size, increment by N
                H5S.get_simple_extent_dims(file_space_id);
                [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
                dims(1) = dims(1)+1; %Appending a single NDArray
                H5D.set_extent (data_id, dims);
                H5S.close(file_space_id);
                start(1) = dims(1)-1;
            end
            H5P.close(lapl_id);
            
            % Get the file space
            file_space_id = H5D.get_space(data_id);
            [~,dims,~] = H5S.get_simple_extent_dims(file_space_id);
            %start(1) = dims(1);
            
            % Select the last N block
            %offset = start;
            H5S.select_hyperslab(file_space_id,'H5S_SELECT_SET',start,ones(rank,1),count,ones(rank,1));
            
            % Mem space
            mem_space_id = H5S.create_simple(rank,count,sizeCount);
           
            
            % Pack the acquisition into the correct struct for writing
            dataOut = data.cplxToStruct();
            
            % Write
            H5D.write(data_id,h5type, ...
                mem_space_id, file_space_id, 'H5P_DEFAULT', dataOut);
            
            % Clean up
            H5S.close(mem_space_id);
            H5S.close(file_space_id);
            H5D.close(data_id);
        end
        

    end
    
end
