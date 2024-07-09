function [dataN] = preproc_behave(rootpath, sname, N)

dataN = cell(N(end),1);
for sbj = N
    if isempty(dir([rootpath sname{sbj} '*.csv']))
        continue
    end
    
    D = dir([rootpath sname{sbj} '*.csv']);
    data = readtable([rootpath D.name]);
    %relo task
    correctness = data.correctness(2:81);
    inputs{1} = data.correctness(2:41); %lo
    inputs{2} = data.correctness(42:81); %re
    answer = data.answer(2:81);
    keys{1} = data.decide_2_keys(find(~cellfun(@isempty,data.decide_2_keys)));
    keys{2} = data.decide_keys(find(~cellfun(@isempty,data.decide_keys)));
    %%%check loss and reward order!!!
    for k = 1:2
        input = inputs{k};
        key = keys{k};
        response = zeros(40,1);
        for i = 1:40
            if input(i) == 1
                if answer{40*(k-1)+i}==key{i}
                    response(i) = 1;
                else
                    response(i)=0;
                end
            elseif input(i) == 0
                if answer{40*(k-1)+i}~=key{i}
                    response(i) = 1;
                else
                    response(i)=0;
                end
            end
        end
        res{k} = response;
    end
    
%     data.input = NaN;
%     data.response = NaN;
    dataN{sbj}.input{1} = inputs{1};
    dataN{sbj}.input{2} = inputs{2};
    dataN{sbj}.response{1} = res{1};
    dataN{sbj}.response{2} = res{2};
end

end
% data 변수에 각 참가자 만큼의 cell 이 만들어지고 각 cell 에 input과 reponse한 0, 1의 값이 담긴 자료가 생성됨
