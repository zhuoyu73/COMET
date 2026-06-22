classdef NDArray < handle
    %NDARRAY Implements ISMRMRD NDArray
    %   Detailed explanation goes here
    
    properties
        data = [];
        dataType = [];
    end
    
    properties(Constant)
        DATA_TYPE = struct( ...
            'USHORT',   uint16(1), ...
            'SHORT',    uint16(2), ...
            'UINT',     uint16(3), ...
            'INT',      uint16(4), ...
            'FLOAT',    uint16(5), ...
            'DOUBLE',   uint16(6), ...
            'CXFLOAT',  uint16(7), ...
            'CXDOUBLE', uint16(8));
    end
    
    
    methods
        
        function obj = NDArray(arg1, dataType)
            switch nargin
                case 0
                    obj.data = [];
                    obj.dataType = obj.DATA_TYPE.CXDOUBLE; %Default
                case 1
                    if(isstruct(arg1))
                        if((isfield(arg1,'real') && isfield(arg1,'imag')))
                            obj.data = complex(arg1.real,arg1.imag);
                        else
                            error('unknown struct input as data to NDArray');
                        end
                    else
                        obj.data = arg1;                 
                    end
                    obj.dataType = obj.getNDArrayDataType(arg1);
                case 2
                    obj.data = arg1;
                    obj.dataType = dataType;
                otherwise
                    error('Incorrect number of arguments supplied to ISMRMRD::NDArray constructor');
            end
        end
        
        function [type,h5type,hltype] = getNDArrayDataType(obj, data)
            if (nargin == 1)
                data = obj.data;
            end
            
            if(isstruct(data))
                if(~(isfield(data,'real') && isfield(data,'imag')))
                    error('unknown structure read in from NDArray');
                end
                if(isa(data.real,'double'));
                    type = obj.DATA_TYPE.CXDOUBLE;
                    h5type = ismrmrd.util.hdf5_datatypes.getType_complexdouble();
                    hltype = 'double';
                end
                
                if(isa(data.real,'single'))
                    type = obj.DATA_TYPE.CXFLOAT;
                    h5type = ismrmrd.util.hdf5_datatypes.getType_complexfloat();
                    hltype = 'single';
                end
            end
            
            if(isa(data,'uint16'));
                type = obj.DATA_TYPE.USHORT;
                h5type = ismrmrd.util.hdf5_datatypes.getType_ushort();
                hltype = 'uint16';
            end
            
            if(isa(data,'int16'));
                type = obj.DATA_TYPE.SHORT;
                h5type = ismrmrd.util.hdf5_datatypes.getType_short();
                hltype = 'int16';
            end
            
            if(isa(data,'uint32'));
                type = obj.DATA_TYPE.UINT;
                h5type = ismrmrd.util.hdf5_datatypes.getType_uint();
                hltype = 'uint32';
            end
            
            if(isa(data,'int32'));
                type = obj.DATA_TYPE.INT;
                h5type = ismrmrd.util.hdf5_datatypes.getType_int();
                hltype = 'int32';
            end
            
            if(isa(data,'double'));
                if(isreal(data))
                    type = obj.DATA_TYPE.DOUBLE;
                    h5type = ismrmrd.util.hdf5_datatypes.getType_double();
                    hltype = 'double';
                else %array is complex
                    type = obj.DATA_TYPE.CXDOUBLE;
                    h5type = ismrmrd.util.hdf5_datatypes.getType_complexdouble();
                    hltype = 'single';
                end
            end
            
            if(isa(data,'single'));
                if(isreal(data))
                    type = obj.DATA_TYPE.FLOAT;
                    h5type = ismrmrd.util.hdf5_datatypes.getType_float();
                    hltype = 'single';
                else %array is complex
                    type = obj.DATA_TYPE.CXFLOAT;
                    h5type = ismrmrd.util.hdf5_datatypes.getType_complexfloat();
                    hltype = 'single';
                end
            end
            
        end
        function strOut = cplxToStruct(obj)
            if ((obj.dataType == obj.DATA_TYPE.CXFLOAT) || (obj.dataType == obj.DATA_TYPE.CXDOUBLE))
                strOut = struct('real',real(obj.data),'imag',imag(obj.data));
            else
                strOut = obj.data;
            end
        end
        
    end
    
end

