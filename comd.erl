-module(comd).
-export([p/0]).

string_to_num(S) ->
	case string:to_float(S) of
    {error,no_float} -> 
    	list_to_integer(S);
     {F,_Rest} -> 
     	io:format("~p~n", [F]),
     	F
    end.

p() ->
    io:format("Starting~n", []),
    Cmd = "./comd\n",
    Port = open_port({spawn,Cmd}, [use_stdio, exit_status]),
    Payload = list_to_binary("1.0 2.0 3.0\n") ,
    io:format("Opened the port: ~w~n", [Port]),
    erlang:port_command(Port, Payload),
    io:format("Sent command to port: ~p~n", [Payload]),
    receive
        {Port, {data, Data}} ->
            io:format("Received data: ~w~n", [Data]),
        	L = string:tokens(Data, " "),
        	io:format("L: ~p~n", [L]),
        	Fun = fun string_to_num/1,
        	lists:foreach(Fun, L),
            io:format("Data: ~p~n", [Data]);
        Other ->
            io:format("Unexpected data: ~p~n", [Other])

    end.