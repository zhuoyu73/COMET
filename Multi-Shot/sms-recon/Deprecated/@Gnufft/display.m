 function display(ob)
%function display(ob)
%	"display" method for this class
Sprintf('Object of class "%s":', class(ob))
disp(struct(ob))
