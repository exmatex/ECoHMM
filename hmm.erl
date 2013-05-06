-module(hmm).

-export([run/1]).

run(N) ->
	{ok, Pid} = comd_sup:start_link(self(), N),
	unlink(Pid),
	AllNames = get_all_names(N),
	run_app(100).

get_all_names(N) when N > 0 ->
	receive
		{comd_started, {Cell, Name}} ->
			io:format("comd_started: [~p]~p~n", [Cell, Name])
	end,
	get_all_names(N-1);
get_all_names(0) ->
	ok.
	
run_app(TS) when TS > 0 ->
	%send_stress(AllNames),
	%get_strain(AllNames),
	%update_stress(AllNames),
	run_app(TS-1);
run_app(0) ->
	ok.
