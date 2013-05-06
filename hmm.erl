-module(hmm).

-export([run/1]).

run(N) ->
	{ok, Pid} = comd_sup:start_link(self()),
	unlink(Pid),
	loop(N).

loop(N) when N > 0 ->
	receive
		{comd_started, {Cell, Name}} ->
			io:format("comd_started: [~p]~p~n", [Cell, Name])
	end,
	loop(N-1);
loop(0) ->
	ok.
	