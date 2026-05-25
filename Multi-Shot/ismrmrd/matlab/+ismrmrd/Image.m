% Image
classdef Image < handle
    
    % Properties
    properties
        
        head = ismrmrd.ImageHeader;
        data = {};
        attribute = {};
        
    end % Properties
    
    % Methods
    methods
        
        function obj = Image(arg1, data, attribute)
            switch nargin
                case 0
                    % No argument constructor
                    % initialize to a single acquisition
                    %extend(obj,1);
                if (length(obj.head) == 1) && isempty(obj.data)

                else
                    error('Object is already initialized! Use other constructor');
                end

                case 1
                    % One argument constructor
                    if ismrmrd.util.isInt(arg1)
                        % First argument is a number
                        M = arg1;
                        extend(obj,M);
                    else
                        % First argument is a header (hopefully)
                        M = length(arg1.version);
                        obj.head = ismrmrd.ImageHeader(arg1);
                        obj.data{M} = [];
                    end
                    
                case 2
                    % Two argument constructor
                    obj.head = ismrmrd.ImageHeader(arg1);
                    M = length(arg1.version);

                    if isempty(data)
                        obj.data{M} = [];
                    else
                        if isreal(data)
                            %dataFromFloat(obj,data); %TODO
                            obj.data{1} = data.data;
                        else
                            obj.data{1} = data.data;
                        end
                    end
                    obj.attribute{1} = '';
                    
                case 3
                    % Three argument constructor
                    obj.head = ismrmrd.ImageHeader(arg1);
                    M = length(arg1.version); %Getting length of header to                    
                    %determine how long the object is.

                    if isempty(data)
                        obj.data{M} = [];
                    else
                        if isreal(data)
                            obj.data{M} = data.data;
                        else
                            obj.data{M} = data.data;
                        end
                    end
                    obj.attribute{M} = attribute; 
                otherwise
                    error('ismrmrd.Image constructor, wrong number of arguments.');
            end
        end
        
        
        function set.head(obj,v)
            obj.head = v;
        end
        
        function set.data(obj,v)
            %obj.data_ = single(complex(v));
            obj.data = v;
        end
        
        function set.attribute(obj,v)
            obj.attribute = v;
        end
        
        function b = isFlagSet(obj,flag)
            bitflag = ismrmrd.FlagBit(flag);
            b = bitflag.isSet(obj.head.flag);
        end
        
        function obj = setFlag(obj,flag)
            bitflag = ismrmrd.FlagBit(flag);
            obj.head.flag = bitor(obj.head.flag, bitflag.bitmask_);
        end
        
        function nacq = getNumber(obj)
            nacq = obj.head.getNumber();
        end
        
        function extend(obj,N)
            % Extend with blank head and empty data.
            if isempty(obj.head)
                M = N;
                obj.head = ismrmrd.ImageHeader(N);
            else
                M = N+obj.getNumber();
                obj.head.extend(N);
            end
            obj.data{M} = [];
            obj.attribute{M} = '';
        end
        
        function append(obj, head, data,attribute)

            Nstart = obj.getNumber + 1;
            Nend   = obj.getNumber + length(head.version);
            Nrange = Nstart:Nend;
            obj.head.append(head);
            if isempty(data) > 0
                obj.data{Nrange} = data;
                if nargin < 4
                    obj.attribute{Nrange} = '';
                else
                    obj.attribute{Nrange} = attribute{Nrange};
                end
            end
            
        end
        
        function img = select(obj, range)
            % Return a copy of a range of acquisitions
            
            % create an empty acquisition
            img = ismrmrd.Image();
            % Fill the header
            img.head = obj.head.select(range);
            % Fill the trajectory and the data
            for p = 1:length(range)
                img.attribute{p} = obj.attribute{range(p)};
                img.data{p} = obj.data{range(p)};
            end
        end
        

    end % Methods
    
end
