 function display(ob)
%function display(ob)
%	"display" method for this class
sprintf('Object of class "%s":', class(ob))
disp(struct(ob))
