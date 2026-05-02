function o = subsref(a, arg)

%     o = builtin('subsref',a,arg);
	if arg.type == '.'
		eval(sprintf('o = a.%s;', arg.subs))
	end
