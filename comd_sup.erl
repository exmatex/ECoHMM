-module(comd_sup).
-behaviour(supervisor).

%% API
-export([start_link/1]).

%% Supervisor callbacks
-export([init/1]).

-define(SERVER, ?MODULE).
-define(COMMAND_LINE, "./comd -x 123.0 -y 234.0 -z 345.0").

%%====================================================================
%% API
%%====================================================================

start_link(Pid) ->
    supervisor:start_link({local, ?SERVER}, ?MODULE, [Pid]).

%%====================================================================
%% Supervisor callbacks
%%====================================================================

init([Pid]) ->
    % Need to dynaically add children with names composed of
    % comd and the cell number of the server.
    Butt =  % Starting one server fo now.
        {tag1,
         {comd_srv, start_link, [comd_0, ?COMMAND_LINE, 0, 0, Pid]},
         transient, 2000, worker, dynamic},
    Nutt =  % Starting one server fo now.
        {tag2,
         {comd_srv, start_link, [comd_1, ?COMMAND_LINE, 1, 0, Pid]},
         transient, 2000, worker, dynamic},
    Futt =  % Starting one server fo now.
        {tag3,
         {comd_srv, start_link, [comd_2, ?COMMAND_LINE, 2, 0, Pid]},
         transient, 2000, worker, dynamic},

    Children = [Butt, Nutt, Futt],
    {ok, {{one_for_one, 3, 10}, Children}}.

%%====================================================================
%% Internal functions
%%====================================================================

